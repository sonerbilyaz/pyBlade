import numpy as np

def Rotate(points, angle, axis):
    
    # Define reversed rotation matrices
    if axis == 'x':
        R = lambda th: np.array([[1.0, 0.0, 0.0],
                                  [0.0, np.cos(th), -np.sin(th)],
                                  [0.0, np.sin(th), np.cos(th)]])
    elif axis == 'y':
        R = lambda th: np.array([[np.cos(th), 0.0, np.sin(th)],
                                  [0.0, 1.0, 0.0],
                                  [-np.sin(th), 0.0, np.cos(th)]])
    elif axis == 'z':
        R = lambda th: np.array([[np.cos(th), -np.sin(th), 0.0],
                                  [np.sin(th), np.cos(th), 0.0],
                                  [0.0, 0.0, 1.0]])
    else:
        raise ValueError("Please enter a valid axis!!")
        
    # Define rotation angles (in radians)
    angle = np.radians(angle)   # Rotation angle
    # Apply rotation
    rotated_points = (points @ R(angle).T)
    
    return rotated_points
