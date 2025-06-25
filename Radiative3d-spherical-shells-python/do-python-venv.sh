#!/bin/bash
# Create a python environment and install the necessary packages.

PYTHON_EXEC="python3"
ENV_NAME="R3Denv"

# Create virtual environment if it doesn't exist
if [ ! -d "$ENV_NAME" ]; then
    echo "Creating virtual environment: $ENV_NAME"
    $PYTHON_EXEC -m venv "$ENV_NAME"
else
    echo "Virtual environment $ENV_NAME already exists."
    echo " Installing the required packages."
    
fi


# Activate environment (for Unix/macOS)
source $ENV_NAME/bin/activate

# Upgrade pip if needed
pip install --upgrade pip

# Install required packages
pip install -r Python/requirements.txt

echo "Environment $ENV_NAME created and dependencies installed."
