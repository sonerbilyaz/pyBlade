import numpy as np

def airfoil_distribution(num_points, dist): 
    if dist == 'cosine_LE-TE':
        # Generate evenly spaced angles between 0 and 90 deg
        angles  = np.linspace(np.radians(0), np.radians(90), num_points//2)
        # Divide the first half (between 0-0.5) into points which cluster towards both ends (0 and 0.5)
        x_half_1, x_half_2= (1-0.5*np.sin(angles)), (0.5*np.sin(angles))

        # Remove repetitive point (midpoint = 0.5)
        x_half_1 = np.delete(np.sort(x_half_1),-1)
        
        # Connect the halves
        points = np.sort(np.hstack((x_half_1,x_half_2)))
        
        # if num_points is NOT an even number
        if num_points/2 != num_points//2:
            Node_ID_one_surf = np.linspace(num_points-2, 1, num_points-2)
            
        # if num_points is an even  number
        if num_points/2 == num_points//2:
            Node_ID_one_surf = np.linspace(num_points-1, 1, num_points-1)
            
    if dist == 'cosine_LE':    
        # Generate evenly spaced angles between 0 and 90 deg
        angles  = np.linspace(np.radians(0), np.radians(90), num_points)
        # Get the distribution
        points = 1-np.sin(angles)
        # Sort in ascending order
        points = np.sort(points)
        
        Node_ID_one_surf = np.linspace(num_points, 1, num_points)
    if dist == 'linear':
        # Remain linear
        points = np.linspace(1, 0, num_points)
        
        Node_ID_one_surf = np.linspace(num_points, 1, num_points)
    
    return points, Node_ID_one_surf.astype(int)

def spanwise_disribution(z_min, z_max, spanwise_panel_num, dist, r_R):
    # Normalizing function when we shift clustering location of the points between (0 and 1)
    def normalize(data, min_val, max_val):
        return min_val + (data - data.min()) * (max_val - min_val) / (data.max() - data.min())
    
    if dist == "cosine":
        # Generate evenly spaced angles between 0 and π
        angles = np.linspace(0, np.pi, spanwise_panel_num)
        # Apply cosine function and normalize to [0, 1]
        points = 0.5 * (1 - np.cos(angles))
        # Map points to [z_min, z_max]
        points = z_min + points * (z_max - z_min)
        return points
    
    if dist == "cosine_TIP":
        angles  = np.linspace(np.radians(0), np.radians(90), int(2.2*(1-r_R)*spanwise_panel_num))
        
        # Uniform part towards root. Remove extra end point
        x1 = np.linspace(0, r_R, int(r_R * spanwise_panel_num))
        x1 = np.delete(x1, -1)
        # cosine TIP part
        x2 = normalize(np.sin(angles), r_R, 1)
        
        # Connect them
        points = np.hstack((x1,x2))
        # Map points to [z_min, z_max]
        points = z_min + points * (z_max - z_min)
        
        ## Remove the 2nd last element since it creates skew panels
        points = np.delete(points,-2)
        
        return points
    
    elif dist == "uniform":
        return np.linspace(z_min, z_max, spanwise_panel_num)      # Linear Distribution
