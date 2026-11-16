#!/usr/bin/env bash
# Run ON Ernie. Default prints plan; --execute starts original-engine tasks.
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
campaign_python="${CAMPAIGN_PYTHON:-/home/s1638308/miniconda3/envs/protac/bin/python3}"
exec "$campaign_python" "$bundle_dir/campaign.py" --root "$bundle_dir" "$@"
