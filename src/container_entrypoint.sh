#!/bin/sh
# Galaxy's Docker runner starts the image as `IMAGE /bin/bash tool_script.sh`,
# so an explicit command runs as given. Flags (or no arguments) keep the
# `docker run IMAGE --point-cloud ...` shorthand for the production runner.
case "${1:-}" in
  "" | -*) exec Rscript /opt/foreststructure/run_julia_memory_safe.R "$@" ;;
  *) exec "$@" ;;
esac
