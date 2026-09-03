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
    if dist == 'uniform':
        # Remain uniform
        points = np.linspace(1, 0, num_points)
        
        Node_ID_one_surf = np.linspace(num_points, 1, num_points)
    
    return points, Node_ID_one_surf.astype(int)

def spanwise_disribution(z_min, z_max, spanwise_panel_num, dist, n=1.2):

    ## First, convert strings from config file into floats
    z_min, z_max, spanwise_panel_num, n = float(z_min), float(z_max), int(spanwise_panel_num), float(n)
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

        ## Remove the 2nd last and 2nd first points since it creates skew panels
        points = np.delete(points,-2)
        points = np.delete(points,[1,2])

        return points
    
    if dist == "polyn_TIP":
        
        ### First, get uniform distribution
        s = np.linspace(0.0, 1.0, spanwise_panel_num)
        ### Then, make it denser towards 1 with power law
        x = 1.0 - (1.0 - s)**n

        # Map points x to [z_min, z_max]
        points = z_min + x * (z_max - z_min)
        
        # ## Remove the last 3rd and 2nd elements since it creates skew panels
        # points = np.delete(points,[-2])
        
        return points
    
    elif dist == "uniform":
        return np.linspace(z_min, z_max, spanwise_panel_num)      # Linear Distribution
