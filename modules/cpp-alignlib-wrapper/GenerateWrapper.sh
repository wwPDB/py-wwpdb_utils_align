#!/bin/bash
#
wrap_pybind11_cli \
    --module_name alignlib \
    --config_path ./alignlib.cfg \
    --output_path ./src \
    --export_json
#