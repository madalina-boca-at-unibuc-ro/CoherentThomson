#!/usr/bin/env bash
# Sync CoherentThomson to a remote build/run server, skipping build artifacts,
# the Python plotting scripts, and local editor/tool caches.
#
# Usage: ./sync_to_remote.sh user@remote:/path/to/CoherentThomson/

set -euo pipefail

DEST="${1:?Usage: $0 user@remote:/path/to/CoherentThomson/}"

rsync -avz \
  --exclude='build/' \
  --exclude='src/build/' \
  --exclude='bin/' \
  --exclude='.venv/' \
  --exclude='.cache/' \
  --exclude='.claude/' \
  --exclude='.vscode/' \
  --exclude='CoherentThomson.code-workspace' \
  /home/madalina/Dropbox/work/bin/c++/codes/CoherentThomson/ "$DEST"
