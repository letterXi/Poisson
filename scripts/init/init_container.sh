#!/bin/bash
echo "START INIT"

sudo chown -R user:user /app

mkdir -p /app/.vscode
mkdir -p /app/.vscode-server

rsync -a --ignore-existing /app/scripts/init/.vscode/*        /app/.vscode
rsync -a --ignore-existing /app/scripts/init/.vscode-server/* /app/.vscode-server

sudo chmod -R 777 /app/.vscode
sudo chmod -R 777 /app/.vscode-server

echo "init_container.sh - DONE"
exec tail -f /dev/null