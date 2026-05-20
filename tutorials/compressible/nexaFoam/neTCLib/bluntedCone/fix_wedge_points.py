#!/usr/bin/env python3
import math
import re
import shutil
import os

# --- Configuration ---
POINTS_FILE = "constant/polyMesh/points"
BACKUP_FILE = "constant/polyMesh/points.bak"

# Based on your logs, the wedge half-angle is 2.0 degrees
WEDGE_HALF_ANGLE = 2.0 

# Tolerance to identify which points belong to the wedge faces (in degrees)
TOLERANCE = 1.0 
# ---------------------

def fix_points():
    if not os.path.exists(POINTS_FILE):
        print(f"Error: {POINTS_FILE} not found. Run this from the case root.")
        return

    # Create a safe backup before modifying
    print(f"Backing up original points file to {BACKUP_FILE}...")
    shutil.copy(POINTS_FILE, BACKUP_FILE)
    
    with open(POINTS_FILE, 'r') as f:
        lines = f.readlines()

    # Regex to capture lines formatted exactly as "(x y z)"
    point_pattern = re.compile(r'^\s*\(\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*\)\s*$')
    
    target_angle_rad = math.radians(WEDGE_HALF_ANGLE)
    
    modified_lines = []
    fixed_count = 0
    
    print("Processing points...")
    for line in lines:
        match = point_pattern.match(line)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))
            
            # Calculate radius in the y-z plane (assuming rotation around x-axis)
            r = math.sqrt(y**2 + z**2)
            
            # Ignore points that are perfectly on the axis of symmetry
            if r > 1e-10: 
                # Current angle in degrees
                theta_deg = math.degrees(math.atan2(z, y))
                
                # Snap to positive wedge angle (wedge1)
                if abs(theta_deg - WEDGE_HALF_ANGLE) < TOLERANCE:
                    y = r * math.cos(target_angle_rad)
                    z = r * math.sin(target_angle_rad)
                    line = f"({x:.16g} {y:.16g} {z:.16g})\n"
                    fixed_count += 1
                    
                # Snap to negative wedge angle (wedge2)
                elif abs(theta_deg - (-WEDGE_HALF_ANGLE)) < TOLERANCE:
                    y = r * math.cos(-target_angle_rad)
                    z = r * math.sin(-target_angle_rad)
                    line = f"({x:.16g} {y:.16g} {z:.16g})\n"
                    fixed_count += 1
            else:
                # Force points on the axis to have exact 0.0 coordinates
                if y != 0.0 or z != 0.0:
                    line = f"({x:.16g} 0 0)\n"
                    fixed_count += 1

        modified_lines.append(line)

    # Write the mathematically perfect coordinates back to the file
    with open(POINTS_FILE, 'w') as f:
        f.writelines(modified_lines)
        
    print(f"Successfully snapped {fixed_count} points to exact mathematical planes.")

if __name__ == "__main__":
    fix_points()
