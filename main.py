#!/usr/bin/env python3
"""photo-geolocator-01 entry: code-only geolocator, no model calls.

Reads the 25 photos of one instance from the env, decodes each JPEG's
DC coefficients (pure stdlib) into a coarse block-colour map, classifies
the scene with hand-written rules into a regional anchor, and posts one
lat/lon answer per photo. Fully deterministic: no clock, no randomness.
"""
import json
import os
import urllib.request

# ---------------------------------------------------------------- JPEG (DC only)

class JpegError(Exception):
    pass


class BitReader:
    """Bits of the entropy segment, with 0xFF00 unstuffing."""

    def __init__(self, data):
        self.d = data
        self.pos = 0
        self.acc = 0
        self.nbits = 0

    def _fill(self):
        while self.nbits <= 24 and self.pos < len(self.d):
            b = self.d[self.pos]
            self.pos += 1
            if b == 0xFF:
                nxt = self.d[self.pos] if self.pos < len(self.d) else 0
                if nxt == 0x00:  # stuffed zero: a real 0xFF byte
                    self.pos += 1
                elif 0xD0 <= nxt <= 0xD7:  # restart marker: realign
                    self.pos += 1
                    self.acc = 0
                    self.nbits = 0
                    return
                else:  # other marker: end of scan data
                    self.pos = len(self.d)
                    return
            self.acc = (self.acc << 8) | b
            self.nbits += 8

    def bits(self, n):
        if n == 0:
            return 0
        self._fill()
        while self.nbits < n:
            self._fill()
            if self.nbits == 0 and self.pos >= len(self.d):
                return 0
        self.nbits -= n
        v = (self.acc >> self.nbits) & ((1 << n) - 1)
        self.acc &= (1 << self.nbits) - 1
        return v

    def align(self):
        self.nbits = 0
        self.acc = 0


class Huff:
    def __init__(self, counts, symbols):
        self.map = {}
        code = 0
        k = 0
        for length in range(1, 17):
            for _ in range(counts[length - 1]):
                self.map[(code, length)] = symbols[k]
                k += 1
                code += 1
            code <<= 1

    def decode(self, br):
        code = 0
        for length in range(1, 17):
            code = (code << 1) | br.bits(1)
            sym = self.map.get((code, length))
            if sym is not None:
                return sym
        raise JpegError("bad huffman code")


def _receive(br, n):
    if n == 0:
        return 0
    v = br.bits(n)
    if v < (1 << (n - 1)):
        v -= (1 << n) - 1
    return v


def decode_dc_grid(data):
    """Return (width, height, rgb_blocks, bw, bh): per-8x8-block mean RGB."""
    if data[:2] != b"\xff\xd8":
        raise JpegError("not a jpeg")
    pos = 2
    quant = {}
    huff = {}
    width = height = 0
    comps = []  # (id, h, v, tq)
    restart = 0
    scan = None
    sof = None
    while pos < len(data):
        if data[pos] != 0xFF:
            pos += 1
            continue
        marker = data[pos + 1]
        if marker == 0xD8:
            pos += 2
            continue
        if marker == 0xD9:
            break
        if 0xD0 <= marker <= 0xD7 or marker == 0x01:
            pos += 2
            continue
        ln = (data[pos + 2] << 8) | data[pos + 3]
        seg = data[pos + 4:pos + 2 + ln]
        if marker == 0xDB:  # DQT
            i = 0
            while i < len(seg):
                pq = seg[i] >> 4
                tq = seg[i] & 15
                i += 1
                table = []
                for _ in range(64):
                    if pq:
                        table.append((seg[i] << 8) | seg[i + 1])
                        i += 2
                    else:
                        table.append(seg[i])
                        i += 1
                quant[tq] = table
        elif marker == 0xC4:  # DHT
            i = 0
            while i < len(seg):
                tc = seg[i] >> 4
                th = seg[i] & 15
                i += 1
                counts = list(seg[i:i + 16])
                i += 16
                total = sum(counts)
                symbols = list(seg[i:i + total])
                i += total
                huff[(tc, th)] = Huff(counts, symbols)
        elif marker in (0xC0, 0xC1, 0xC2):  # SOF
            sof = marker
            height = (seg[1] << 8) | seg[2]
            width = (seg[3] << 8) | seg[4]
            ncomp = seg[5]
            comps = []
            for c in range(ncomp):
                o = 6 + c * 3
                comps.append((seg[o], seg[o + 1] >> 4, seg[o + 1] & 15, seg[o + 2]))
        elif marker == 0xDD:  # DRI
            restart = (seg[0] << 8) | seg[1]
        elif marker == 0xDA:  # SOS: scan data follows
            ns = seg[0]
            sel = []
            for c in range(ns):
                o = 1 + c * 2
                sel.append((seg[o], seg[o + 1] >> 4, seg[o + 1] & 15))
            scan = data[pos + 2 + ln:]
            break
        pos += 2 + ln
    if sof == 0xC2:
        raise JpegError("progressive jpeg not supported")
    if not comps or scan is None:
        raise JpegError("missing frame or scan")
    hmax = max(c[1] for c in comps)
    vmax = max(c[2] for c in comps)
    mcux = (width + 8 * hmax - 1) // (8 * hmax)
    mcuy = (height + 8 * vmax - 1) // (8 * vmax)
    br = BitReader(scan)
    # per-component block grids of dequantised DC averages
    grids = []
    for (cid, h, v, tq) in comps:
        bw = (width * h + 8 * hmax - 1) // (8 * hmax)
        bh = (height * v + 8 * vmax - 1) // (8 * vmax)
        grids.append([[0.0] * bw for _ in range(bh)])
    dc_pred = [0] * len(comps)
    comp_by_id = {c[0]: i for i, c in enumerate(comps)}
    sel_idx = [(comp_by_id[cid], huff[(0, tdc)], huff[(1, tac)],
                quant[comps[comp_by_id[cid]][3]][0])
               for (cid, tdc, tac) in sel]
    total_mc = mcux * mcuy
    for m in range(total_mc):
        if restart and m > 0 and m % restart == 0:
            br.align()
            dc_pred = [0] * len(comps)
        my, mx = divmod(m, mcux)
        for (ci, htab, atab, q0) in sel_idx:
            h, v = comps[ci][1], comps[ci][2]
            for by in range(v):
                for bx in range(h):
                    s = htab.decode(br)
                    dc_pred[ci] += _receive(br, s)
                    avg = dc_pred[ci] * q0 / 8.0 + 128.0
                    gx = mx * h + bx
                    gy = my * v + by
                    if gy < len(grids[ci]) and gx < len(grids[ci][0]):
                        grids[ci][gy][gx] = avg
                    # skip AC coefficients
                    n = 1
                    while n < 64:
                        s = atab.decode(br)
                        if s == 0x00:
                            break
                        run = s >> 4
                        size = s & 15
                        if s == 0xF0:
                            n += 16
                            continue
                        n += run + 1
                        if size:
                            br.bits(size)
    # upsample chroma grids to the luma grid, convert YCbCr -> RGB
    lw = len(grids[0][0])
    lh = len(grids[0])
    rgb = [[(0, 0, 0)] * lw for _ in range(lh)]
    if len(comps) == 1:
        for y in range(lh):
            for x in range(lw):
                g = max(0, min(255, int(grids[0][y][x])))
                rgb[y][x] = (g, g, g)
    else:
        for y in range(lh):
            for x in range(lw):
                yy = grids[0][y][x]
                cb = grids[1][min(y * len(grids[1]) // lh, len(grids[1]) - 1)][
                    min(x * len(grids[1][0]) // lw, len(grids[1][0]) - 1)] - 128.0
                cr = grids[2][min(y * len(grids[2]) // lh, len(grids[2]) - 1)][
                    min(x * len(grids[2][0]) // lw, len(grids[2][0]) - 1)] - 128.0
                r = max(0, min(255, int(yy + 1.402 * cr)))
                g = max(0, min(255, int(yy - 0.34414 * cb - 0.71414 * cr)))
                b = max(0, min(255, int(yy + 1.772 * cb)))
                rgb[y][x] = (r, g, b)
    return width, height, rgb, lw, lh


# ---------------------------------------------------------------- features

def _frac(sel, pred):
    if not sel:
        return 0.0
    return sum(1 for p in sel if pred(p)) / len(sel)


def _band(rgb, bw, bh, y0, y1, x0=0.0, x1=1.0):
    out = []
    for y in range(int(bh * y0), int(bh * y1)):
        for x in range(int(bw * x0), int(bw * x1)):
            out.append(rgb[y][x])
    return out


def features(rgb, bw, bh):
    flat = [p for row in rgb for p in row]
    n = len(flat)
    green = lambda p: p[1] > p[0] + 12 and p[1] > p[2] + 12
    blue = lambda p: p[2] > p[0] + 15 and p[2] > p[1] + 8
    tan = lambda p: p[0] > p[2] + 25 and p[0] > 110 and p[1] > 70
    reddirt = lambda p: p[0] > 110 and p[0] > p[1] + 30 and p[1] > p[2] + 10
    bright = lambda p: p[0] > 200 and p[1] > 200 and p[2] > 200
    gray = lambda p: abs(p[0] - p[1]) < 12 and abs(p[1] - p[2]) < 12
    return {
        "green": _frac(flat, green),
        "blue": _frac(flat, blue),
        "tan": _frac(flat, tan),
        "reddirt": _frac(flat, reddirt),
        "bright": _frac(flat, bright),
        "gray": _frac(flat, gray),
        "sky_blue": _frac(_band(rgb, bw, bh, 0.0, 0.30), blue),
        "sky_bright": _frac(_band(rgb, bw, bh, 0.0, 0.30), bright),
        "bot_tan": _frac(_band(rgb, bw, bh, 0.55, 1.0), tan),
        "bot_red": _frac(_band(rgb, bw, bh, 0.55, 1.0), reddirt),
        "bot_green": _frac(_band(rgb, bw, bh, 0.55, 1.0), green),
        "mid_green": _frac(_band(rgb, bw, bh, 0.30, 0.70), green),
        "lum": sum(p[0] + p[1] + p[2] for p in flat) / n / 3,
    }


# ---------------------------------------------------------------- classifier

# Anchors are round representative points for a region, general knowledge.
PRIOR = (56.0, 10.0)


def classify(f):
    arid = f["tan"] > 0.30
    lush = f["green"] > 0.22
    red = f["bot_red"] > 0.05
    snow = f["bright"] > 0.25
    skyblue = f["sky_blue"] > 0.45
    if arid and f["green"] < 0.05:
        return (-24.0, -66.0)  # arid Andes / dry mountains
    if arid:
        return (-15.0, -50.0)  # dry interior, Brazil-like
    if snow and f["green"] < 0.08:
        return (58.0, 30.0)  # snow / bright winter
    if lush:
        if red:
            return (-5.0, 35.0)  # lush + red earth: equatorial Africa
        return (0.0, 105.0)  # lush tropics: SE Asia / Oceania
    if red:
        return (-8.0, 30.0)  # red dirt tropics
    if skyblue and f["green"] < 0.08:
        return (60.0, 25.0)  # clear northern skies
    if skyblue:
        return (52.0, 10.0)  # temperate, blue sky
    if f["green"] > 0.08:
        return (48.0, 8.0)  # temperate green
    return PRIOR


# ---------------------------------------------------------------- env client

def http_json(url, data=None):
    req = urllib.request.Request(url, data=data)
    if data is not None:
        req.add_header("Content-Type", "application/json")
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def http_bytes(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def main():
    env_url = os.environ.get("MAKEFUN_ENV_URL", "").rstrip("/")
    if not env_url:
        raise SystemExit("MAKEFUN_ENV_URL is not set")
    listing = http_json(env_url + "/photos")
    answers = []
    for photo in listing["photos"]:
        pid = photo["id"]
        lat, lon = PRIOR
        try:
            data = http_bytes(env_url + photo["url"])
            _w, _h, rgb, bw, bh = decode_dc_grid(data)
            lat, lon = classify(features(rgb, bw, bh))
        except Exception:
            lat, lon = PRIOR
        lat = max(-90.0, min(90.0, lat))
        lon = max(-180.0, min(180.0, lon))
        answers.append({"id": pid, "lat": round(lat, 6), "lon": round(lon, 6)})
    http_json(env_url + "/answers",
              data=json.dumps({"answers": answers}).encode("utf-8"))


if __name__ == "__main__":
    main()
