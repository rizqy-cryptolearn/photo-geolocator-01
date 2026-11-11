# photo-geolocator-01 — code-only DC-colour geolocator

A model-free entry for make.fun's photo-geolocator-01: it places each
street-level photo by reading the photo's own pixels — nothing else.

## How it works

1. Fetches the instance's 25 photos from the env (`GET /photos`, then one
   `GET /photos/<id>.jpg` each).
2. Decodes every JPEG's **DC coefficients only** (pure Python standard
   library, no image libraries in the sandbox) into a coarse per-8x8-block
   colour map.
3. Measures scene fractions from that map — vegetation, sky blue, tan dry
   ground, red earth, snow-bright, grey asphalt — top/middle/bottom bands.
4. Applies hand-written rules that map a scene to a regional anchor
   (round representative coordinates: arid Andes, lush tropics, snow
   winter, temperate green, and a global prior). No model calls, no
   network beyond the env, no lookup of any photo's place.
5. Posts all 25 answers in one `POST /answers` and exits.

Deterministic: no clock, no randomness, every answer is a function of the
photo bytes alone, so the live run and the replay agree exactly.

## Run locally

The kit's env server serves a practice instance on a port; point the agent
at it:

```sh
MAKEFUN_ENV_DATA_DIR=<kit>/practice/trial-0/input/env \
MAKEFUN_STATE_DIR=/tmp/geo-state MAKEFUN_ENV_PORT=18080 \
node <kit>/program/env/server.mjs &
MAKEFUN_ENV_URL=http://127.0.0.1:18080 python3 main.py
```
