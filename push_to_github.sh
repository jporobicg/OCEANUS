#!/bin/bash

# Simple script to connect, pull, push, and disconnect from GitHub
cd /home/por07g/Documents/Code_Tools/OCEANUS

echo "1) Connecting to GitHub..."
git remote add origin git@github.com:jporobicg/OCEANUS.git 2>/dev/null || git remote set-url origin git@github.com:jporobicg/OCEANUS.git

echo "2) Pulling latest changes from GitHub..."
git pull origin main

echo "3) Pushing to GitHub..."
git push -u origin main

echo "4) Disconnecting from GitHub..."
git remote remove origin

echo "Done! Repository pulled, pushed, and disconnected."
