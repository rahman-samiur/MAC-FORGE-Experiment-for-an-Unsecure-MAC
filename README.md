# MAC-FORGE-Experiment-for-an-Unsecure-MAC
This project was a requirement for the CS549 (Cryptography) course

# Python Project

## Overview

This project consists of a Python application with the main entry point being `main.py`, located in the source code folder along with several other supporting files.

## Requirements

- Python 3.x
- No additional external libraries required beyond the Python standard library
- The project only uses built-in modules:
  - `random`
  - `os`
  - `json`

## Project Structure

```
project/
│
├── main.py          # Main entry point for the application
└── [other files]    # Supporting files used by the application
```

## Execution Environment

This application is designed to run in any environment with Python 3.x installed. No virtual environment or package installation is required since it only uses standard library modules.

## Running the Application

To run the application, navigate to the project directory and execute:

```bash
python main.py
```

If you're on a system where Python 3 is not the default Python version, use:

```bash
python3 main.py
```

## Configuration

Any configuration settings are likely stored in JSON format and loaded by the application. Check the source files for more information about specific configuration options.

## Troubleshooting

- Ensure Python 3.x is installed on your system
- Verify that the project files have appropriate read/write permissions
- If the application attempts to write files, make sure you have write permissions in the directory
