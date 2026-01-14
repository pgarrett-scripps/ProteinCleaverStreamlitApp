# List available commands
default:
    @just --list

# Run the streamlit app
run:
    uv run streamlit run app.py

# Install dependencies
install:
    uv venv
    uv pip install -r requirements.txt

sync:
    uv pip install -r requirements.txt

# Upgrade requirements.txt to newest versions using uv
upgrade-deps:
    @echo "Updating requirements.txt..."
    uv pip compile requirements.txt -o requirements.txt --upgrade

# Install dev dependencies
install-dev:
    uv pip install -r requirements-dev.txt

# Run code formatting
format:
    uv run ruff format .

# Run linter
lint:
    uv run ruff check .

# Run type checking
ty:
    uv run ty check .
