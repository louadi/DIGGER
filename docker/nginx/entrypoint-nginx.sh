#!/usr/bin/env bash

# ======================================================================================
# title			:  entrypoint-nginx.sh
# description	:  This script replaces specific variables through environmental variables
# usage			:  bash entrypoint-nginx.sh
# author       	:  Kevin Yuan
# date			:  24.01.2020
# version		:  1.0
# notes			:  This script is called by the docker file
# ======================================================================================

set -eu

envsubst '${NGINX_PUBLISHED_PATH}' < /etc/nginx/conf.d/nginx.conf.template > /etc/nginx/conf.d/nginx.conf

exec "$@"