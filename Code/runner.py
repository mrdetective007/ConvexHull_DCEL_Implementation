"""
This file is used to run the generator.py, main.cpp, and visualization.py files.
"""
#!/opt/homebrew/bin/python3
import os

# Run generator.py
os.system('python3 generator.py > input.txt')

# compile main.cpp
os.system('g++-12 -std=c++17 main.cpp -o main')

# Run main.cpp
os.system('./main')

# Run visualizer.py
os.system('python3 visualization.py')
