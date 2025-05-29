#!/bin/bash

# Get total system memory in MB
total_mem=$(free -m | awk '/^Mem:/{print $2}')
available_mem=$(free -m | awk '/^Mem:/{print $7}')

echo "System Memory Information:"
echo "Total Memory: ${total_mem}MB"
echo "Available Memory: ${available_mem}MB"
echo "Required Memory (approximate): 70MB"

# Check if we have enough memory
if [ $available_mem -lt 70 ]; then
    echo "WARNING: Less than 70MB of memory available. This might cause issues."
fi

# Get stack size limit
stack_size=$(ulimit -s)
echo "Stack Size Limit: ${stack_size}KB"

# Check if stack size might be too small
if [ $stack_size -lt 8192 ]; then
    echo "WARNING: Stack size is less than 8MB. Consider increasing it with:"
    echo "ulimit -s 8192"
fi 