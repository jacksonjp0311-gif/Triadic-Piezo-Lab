#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
if command -v python3 >/dev/null 2>&1; then
  python3 -m http.server 8765 --bind 127.0.0.1 >/tmp/triad-lab-server.log 2>&1 &
  sleep 1
  URL="http://127.0.0.1:8765/OPEN_LAB.html"
else
  URL="file://$(pwd)/OPEN_LAB.html"
fi
if command -v open >/dev/null 2>&1; then open "$URL"; elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL"; else printf '%s\n' "$URL"; fi
