#!/bin/bash
export HISTIGNORE="*passwd*"
echo "${1}:${2}" | sudo chpasswd