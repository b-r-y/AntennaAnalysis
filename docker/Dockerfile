# Use the latest Alpine Linux image
FROM alpine:latest

# Update and install dependencies
RUN apk update && apk add --no-cache python3 py3-pip git texmf-dist texlive-full

# Set up a working directory
WORKDIR /workspace

# Create and activate a Python virtual environment
RUN python3 -m venv /workspace/.venv

# Set the PATH for using the virtual environment
ENV PATH="/workspace/.venv/bin:$PATH"

# Copy requirements.txt (assuming it exists in workspace)
COPY requirements.txt /workspace/

# Install Python dependencies inside the virtual environment
RUN pip install --upgrade pip && pip install -r requirements.txt

# Default command (optional)
CMD ["sh"]
