import numpy as np
import xarray as xr

class DinoConfiguration:
    def __init__(
            self,
            #grid
            res     = 1.0,              # Horizontal resolution at equator [degrees]
            r_earth = 6371229.,         # radius of the earth
            lon_min = 0.,               # Minimum longitude on inner U-point
            lon_max = 50.,              # Maximum longitude on inner U-point
            lat_min = -70.,             # Minimum latitude on V-point (approx.)
            lat_max = 70.,              # Maximum latitude on V-point (approx.)
            #bathymet   
            lat_channel_min = -65.,     # Minimum channel latitude on T-point (approx.)     
            lat_channel_max = -45.,     # Maximum channel latitude on T-point (approx.)
            slope = 3.,                 # slope
            H_max = 4000.,              # Maximum depth of the bathymetry on W-point
            H_min = 2000.,              # Minimum depth of the bathymetry on W-point (approx.)
            #vertical coordinate    
            K = 36,                     # Number of vertical levels
            dz_min = 10.,               # Minimum value of e3 at the surface [meters]
            k_th = 35.,                 # Position of the inflexion point
            h_co = 1000.,               # Depth of the connection between z- and s-coordinates [meters]
            a_cr = 10.5,                # Slope of the tanh function
            hybrid_coord = False,       # Switch to hybrid sigma-z-coordinate
            #drake_si   
            slope_ds = 4.,              # Length scale of the slope of the sill [degrees]
            depth_ds = 2500.,           # Depth of the sill [meters]
            #forcing
            utau_nodes = [              # nodes of latitude for windstress profile
                -45. ,-15. ,0. ,15. ,45.   
            ],
            utau_values = [             # values of windstress profile at utau_nodes
                0.2 ,-0.1 ,-0.02 ,-0.1 ,0.1
            ],
            t_star_n = 5.0,             # temperature restoring at northern boundary
            t_star_s = -0.5,            # temperature restoring at southern boundary
            t_star_eq = 27.,            # temperature restoring at equator
            s_star_n = 35.,             # salinity restoring at northern boundary
            s_star_s = 35.1,            # salinity restoring at southern boundary
            s_star_eq = 37.25,          # salinity restoring at equator


    ):
        self.res                = res
        self.r_earth            = r_earth
        self.lon_min            = lon_min
        self.lon_max            = lon_max
        self.lat_min            = lat_min
        self.lat_max            = lat_max
        self.lat_channel_min    = lat_channel_min
        self.lat_channel_max    = lat_channel_max
        self.slope              = slope
        self.H_max              = H_max
        self.H_min              = H_min
        self.K                  = K
        self.dz_min             = dz_min
        self.k_th               = k_th
        self.h_co               = h_co
        self.a_cr               = a_cr
        self.hybrid_coord       = hybrid_coord
        self.slope_ds           = slope_ds
        self.depth_ds           = depth_ds
        self.utau_nodes         = utau_nodes
        self.utau_values        = utau_values
        self.t_star_n           = t_star_n 
        self.t_star_s           = t_star_s 
        self.t_star_eq          = t_star_eq
        self.s_star_n           = s_star_n 
        self.s_star_s           = s_star_s 
        self.s_star_eq          = s_star_eq

        # set coordinates:
        self.lon_u, self.lon_t  = self.get_longitude()
        self.lat_v, self.lat_t  = self.get_latitude()

        # set meshgrid
        self.lon_mesh_t, self.lat_mesh_t = np.meshgrid(self.lon_t, self.lat_t)
        self.lon_mesh_f, self.lat_mesh_f = np.meshgrid(self.lon_u, self.lat_v)

        # set bathymetry:
        self.bathymetry         = self.get_bathymetry()

        # add sill in drake-passage
        self.bathymetry = self.add_gauss_ring() 

        # set vertical coordinate
        self.depth_w, self.depth_t = self.get_vertical_coordinate()

        # set scale factors
        self.e1t, self.e1f, self.e3t, self.e3w = self.get_scale_factors()

        # rest of scale factors are identical to previous
        self.e2t = self.e1t.copy()
        self.e2f = self.e1f.copy()
        self.e1u = self.e1f.copy()
        self.e2u = self.e2t.copy()
        self.e1v = self.e1t.copy()
        self.e2v = self.e2f.copy()

        if self.hybrid_coord: # TODO: implement vertical coordinates on u-, v-, f-points
            print("Hybrid vertical coordinate selected, but adapted e3u/v/f is not implemented yet")
            self.e3u = None
            self.e3v = None
            self.e3t = None
        else:
            self.e3u = self.e3t.copy()
            self.e3v = self.e3t.copy()
            self.e3f = self.e3t.copy()
        
        #set surface boundary conditions
        self.utau       = self.get_utau()
        self.q_solar    = self.get_q_solar()
        self.t_star     = self.get_T_star()
        self.s_star     = self.get_S_star()
        self.rho_star   = self.get_rho_star()

    def get_longitude(self):
        '''
        Compute the longitude coordinate.

        Returns:
        - lon_u:    longitude coordinate on u-points
        - lat_t:    longitude coordinate on t-points
        '''
        lon_u = np.arange(self.lon_min, self.lon_max + 2 * self.res, self.res)
        lon_t = np.arange(self.lon_min - self.res / 2, self.lon_max + self.res, self.res)
        return(lon_u, lon_t)

    def get_latitude(self):
        '''
        Compute the latitude coordinate.

        Returns:
        - lat_v:    latitude coordinate on u-points
        - lat_t:    latitude coordinate on t-points
        '''
        rad     = ( np.pi / 180 )
        j_eq_n  = self.mercator_projection(self.lat_max, self.res)
        j_eq_s  = self.mercator_projection(self.lat_min, self.res)

        N_y     = (j_eq_s - j_eq_n) + 1
        j_t = np.arange(j_eq_n, j_eq_s + 1, 1.0)
        j_v = j_t + 0.5
        lat_v   = np.arcsin(np.tanh(self.res * rad * j_v)) / rad
        lat_t   = np.arcsin(np.tanh(self.res * rad * j_t)) / rad
        return(lat_v, lat_t)
    
    def get_scale_factors(self):
        """
        Compute the scale factors e1t, e1f, e3t, e3w.
        These are enough to define the scale factors on all grid-points.

        Returns:
        - e1t, e1f, e3t, e3w
        """
        rad = ( np.pi / 180.)

        # e1
        e1t = self.r_earth * rad * np.cos( rad * self.lat_mesh_t ) * self.res
        e1f = self.r_earth * rad * np.cos( rad * self.lat_mesh_f ) * self.res
        # e3w
        e3w = 2. * ( self.depth_t - self.depth_w ) 
        e3w[:,:,1:] = self.depth_t[:,:,1:] - self.depth_t[:,:,:-1]
        # e3t
        e3t = 2. * ( self.depth_t - self.depth_w ) 
        e3t[:,:,:-1] = self.depth_w[:,:,1:] - self.depth_w[:,:,:-1]

        return(e1t, e1f, e3t, e3w)
    
    def get_mask(self):
        pass

        
    def get_vertical_coordinate(self):
        """
        Define vertical s-coordinate system using the 2D bathymetry field.

        Returns:
        - depth_w:  Depth of the levels at w points 
                    (3D numpy array with shape (bathy.shape[0], bathy.shape[1], K))
        - depth_t:  Depth of the levels at t points 
                    (3D numpy array with shape (bathy.shape[0], bathy.shape[1], K))
        """
        if self.hybrid_coord:           # vertical coordinate follows bathymetry below h_co
            bottom = self.bathymetry
        else:                           # flat bottom
            bottom = np.ones_like(self.bathymetry) * self.H_max

        # 1D profile
        depth_1d = self.mi96_1d(self.H_max, self.dz_min, 0.0, self.k_th, 0, self.a_cr, self.K)
        depthw_1d = depth_1d[0, :]
        deptht_1d = depth_1d[1, :]

        # Nearest index and its depth value to the depth where s-coordinates should start
        k_const = np.argmin(np.abs(depthw_1d - self.h_co))

        # compute e3w from depth:
        e3w_1d     = 2. * ( deptht_1d - depthw_1d ) 
        e3w_1d[1:] = deptht_1d[1:] - deptht_1d[:-1]

        depth_top = np.broadcast_to(depth_1d[:,:k_const], (self.bathymetry.shape[0], self.bathymetry.shape[1], 2, k_const))
        
        K_top = self.K - k_const

        # vectorize mi96_1d
        vec_mi96_1d = np.vectorize(
            lambda bathy_val: self.mi96_1d(bathy_val, e3w_1d[k_const], depthw_1d[k_const], self.k_th, k_const, self.a_cr, self.K),
            signature='()->(2,K_top)'
        )

        # compute the 3D depth levels at t and w points
        depth_bottom  = vec_mi96_1d(bottom)
        depth = np.concatenate([depth_top, depth_bottom], axis=3)
        depth_w       = depth[:,:,0,:]
        depth_t       = depth[:,:,1,:]
        return(depth_w, depth_t)

    def add_gauss_ring(self):
        """
        Add a ring of Gaussian bumps to the bathymetry representing an idealized drake-sill.

        Returns:
        float: Resulting updated bathymetry [meters]
        """
        lat_0       = (self.lat_channel_max + self.lat_channel_min) / 2
        lon_0       = self.lon_t.min()
        radius      = abs(self.lat_channel_max - self.lat_channel_min) / 2

        # taper the gaussian ring eastward 
        taper       = self.smooth_step(self.lon_mesh_t, lon_0, lon_0 + self.slope_ds)
        # shift grid
        zx = self.lon_mesh_t - lon_0
        zy = self.lat_mesh_t - lat_0

        # Calculate the exponent part of the Gaussian function
        exp_arg = (- zx**2 - zy**2 + 2 * radius * np.sqrt(zx**2 + zy**2) - radius**2) / self.slope_ds**2

        # Calculate the resulting bathymetry depth
        gauss_ring = np.where(self.bathymetry >= self.depth_ds,
                         (self.depth_ds - self.bathymetry) * np.exp(exp_arg) + self.bathymetry,
                         self.bathymetry)
        
        # update bathymetry
        new_bathymetry = taper * gauss_ring + (1. - taper) * self.bathymetry
        return(new_bathymetry)

    def get_bathymetry(self):
        '''
        Compute bathymetry as 2D numpy array.
        '''
        width           = abs(self.lon_max - self.lon_min)
        channel_width   = abs(self.lat_channel_max - self.lat_channel_min)

        zy_cha = self.exp_bathymetry(self.lat_mesh_t, self.lat_channel_min, self.lat_channel_max, width, self.slope, channel_width / 2)

        zx_raw = self.exp_bathymetry(self.lon_mesh_t, self.lon_min, self.lon_max, width, self.slope, channel_width / 2)
        zx = zx_raw * (1.0 - zy_cha) + zy_cha

        slope_lat = np.cos( np.pi * self.lat_max /180) * self.slope
        zy = self.exp_bathymetry(self.lat_mesh_t, self.lat_min, self.lat_max, width, slope_lat, channel_width / 2)

        pbathy = zx * zy * (self.H_max - self.H_min) + self.H_min

        return pbathy

    def get_utau(self):
        '''
        Compute zonal wind stress as 1D meridional profile.
        '''
        #Append/prepend zero wind-stress at the boundary
        utau_nodes  = [self.lat_t[0], *self.utau_nodes, self.lat_t[-1]]
        utau_values = [0., *self.utau_values, 0.]
        utau = self.tau_from_nodes(utau_nodes, utau_values, self.lat_t)
        return(utau)
    
    def get_q_solar(self, day_of_year=0.0):
        '''
        Compute solar heatflux as meridional profile.
        '''
        q_solar = np.maximum(
            230. * np.cos(
                np.pi * (self.lat_t - 23.5 * np.expand_dims(self.seasonal_cycle(day_of_year, lag=0.), axis=-1)) / 180.
            ),
            0.
        )
        return(q_solar)
    
    def get_S_star(self):
        '''
        Compute salinity restoring as 1D meridional profile.
        '''
        # Calculate common terms
        cos_term = np.cos(2 * np.pi * self.lat_t / (self.lat_max - self.lat_min))
        exp_term = np.exp(-self.lat_t ** 2 / 7.5 ** 2)

        # Calculate S_star depending on hemisphere
        s_star = np.where(self.lat_t <= 0,
                          self.s_star_s + (self.s_star_eq - self.s_star_s) * (1 + cos_term) / 2 - 1.25 * exp_term,
                          self.s_star_n + (self.s_star_eq - self.s_star_n) * (1 + cos_term) / 2 - 1.25 * exp_term)
        return s_star
    
    def get_T_star(self, day_of_year=0.0):
        '''
        Compute salinity restoring as meridional profile.
        '''
        t_star_n = self.t_star_n + 3.  * self.seasonal_cycle(day_of_year, lag=30.)
        t_star_s = self.t_star_s + 0.5 * self.seasonal_cycle(day_of_year, lag=30.)

        #in case day_of_year is an array for broadcasting
        t_star_s = np.expand_dims(t_star_s, axis=-1)  # Shape becomes (360, 1)
        t_star_n = np.expand_dims(t_star_n, axis=-1)  # Shape becomes (360, 1)
        
        # Calculate common term
        sin_term = np.sin(np.pi * (self.lat_t + self.lat_max) / (self.lat_max - self.lat_min))

        # Calculate S_star depending on hemisphere
        t_star = np.where(
            self.lat_t <= 0,
            t_star_s + (self.t_star_eq - t_star_s) * sin_term,    # Southern hemisphere
            t_star_n + (self.t_star_eq - t_star_n) * sin_term     # Northern hemisphere
        )
        return t_star
    
    def get_rho_star(self, day_of_year=0.0):
        '''
        Compute density restoring as meridional profile.

        Use simplified equation of state.
        '''
        s_star = self.s_star
        t_star = self.get_T_star(day_of_year)

        rho_star = (
                    - 0.165 * (1. + 0.5 * 0.06 * ( t_star - 10.)) * ( t_star - 10.)
                    + 0.7655 * ( s_star - 35.)
                ) + 1026
        return rho_star

    def to_xarray(self):
        '''
        Write all configuration fields to an xarray dataset.
        '''
        ds = xr.Dataset(
            data_vars=dict(
                e1t         = (['y_c', 'x_c'], self.e1t),
                e1u         = (['y_c', 'x_f'], self.e1u),
                e1v         = (['y_f', 'x_c'], self.e1v),
                e1f         = (['y_f', 'x_f'], self.e1f),
                e2t         = (['y_c', 'x_c'], self.e2t),
                e2u         = (['y_c', 'x_f'], self.e2u),
                e2v         = (['y_f', 'x_c'], self.e2v),
                e2f         = (['y_f', 'x_f'], self.e2f),
                e3t         = (['y_c', 'x_c', 'z_c'], self.e3t),
                e3u         = (['y_c', 'x_f', 'z_c'], self.e3u),
                e3v         = (['y_f', 'x_c', 'z_c'], self.e3v),
                e3f         = (['y_f', 'x_f', 'z_c'], self.e3f),
                e3w         = (['y_c', 'x_c', 'z_f'], self.e3f),
                bathymetry  = (['y_c', 'x_c'], self.bathymetry),
                depth_t     = (['y_c', 'x_c', 'z_c'], self.depth_t),
                depth_w     = (['y_c', 'x_c', 'z_f'], self.depth_w),
                utau        = (['y_c'], self.utau),
                q_solar     = (['y_c'], self.q_solar),
                t_star      = (['y_c'], self.t_star),
                s_star      = (['y_c'], self.s_star),
                rho_star    = (['y_c'], self.rho_star),
            ),
            coords=dict(
                lon_c       = ('x_c', self.lon_t),
                lon_f       = ('x_f', self.lon_u),
                lat_c       = ('y_c', self.lat_t),
                lat_f       = ('y_f', self.lat_v),
                depth_t_1d  = ('z_c', self.depth_t[round(len(self.lat_t) / 2), round(len(self.lon_t) / 2), :]),
                depth_w_1d  = ('z_f', self.depth_w[round(len(self.lat_t) / 2), round(len(self.lon_t) / 2), :]),
            ),
            attrs=dict(
                description='DINO configuration created with DinoConfiguration class',
                resolution=self.res
            )
        )
        return(ds)
    
    @staticmethod
    def mercator_projection(lat, res):
        '''
        Compute the number of gridpoints to the equator from a given latitude (Mercator projection).

        Parameters:
        - phi: approximate latitude in degrees
        - res: resolution of the grid in degree

        Returns:
        - N: number of gridpoints between equator and approximate latitude
        '''
        arg = np.pi / 4.0 - np.pi / 180.0 * lat / 2.0
        jeq = np.abs(180.0 / np.pi * np.log(np.cos(arg) / np.sin(arg)) / res)

        if lat > 0:
            jeq = -jeq

        N = np.round(jeq).astype(int)
        return(N)

    @staticmethod
    def mi96_1d(bathy, dz_min, h_co, k_th, k_const, a_cr, K):
        '''
        Calculate depth levels for a given set of input parameters using Madec & Imbard (1996) function.

        Parameters:
        - bathy:    Full depth of the water column
        - dz_min:   Value of e3 at the connection point [meters]
        - h_co:     Depth of the connection between z- and s-coordinates [meters]
        - k_th:     Position of the inflexion point
        - k_const:  Number of levels with z-coordinates [meters]
        - a_cr:     Slope of the tanh function
        - K:        Number of vertical levels

        Returns:
        - depth:    Depth of the levels (2D numpy array with shape (2, K - k_const))
        '''

        depth = np.zeros((2, K - k_const))

        # Compute parameters za0, za1, and zsur
        a0 = (dz_min - (bathy - h_co) / (K - 1 - k_const)) / (
            np.tanh((1 - k_th) / a_cr) -
            a_cr / (K - 1 - k_const) * (np.log(np.cosh((K - k_const - k_th) / a_cr)) -
                                      np.log(np.cosh((1 - k_th) / a_cr)))
        )

        a1 = dz_min - a0 * np.tanh((1 - k_th) / a_cr)
        a2 = -a1 - a0 * a_cr * np.log(np.cosh((1 - k_th) / a_cr))

        # Calculate depth at T and W-points
        for jk in range(1, K - k_const + 1):
            zw = float(jk)
            zt = float(jk) + 0.5

            depth[0, jk - 1] = a2 + a1 * zw + a0 * a_cr * np.log(np.cosh((zw - k_th) / a_cr)) + h_co
            depth[1, jk - 1] = a2 + a1 * zt + a0 * a_cr * np.log(np.cosh((zt - k_th) / a_cr)) + h_co

        return depth

    @staticmethod
    def smooth_step( lam, lam_left, lam_right):
        '''
        Compute a smooth step function for tapering

                    |                                ***|*****---> 1 
                    |                           *****   | 
                    |                       ****        | 
                    |                     **            | 
                    |                    *              | 
                    |                   *               | 
                    |                 **                | 
                    |              ***                  | 
                    |         *****                     | 
                    |   ******                          |
         0 <---*****|***                                |
                    v                                   v
                 lam_left                            lam_right

        Parameters:
        - lam:          coordinate value (e.g. lon or lat)
        - lam_left:     lower end of smoothing function
        - lam_right:    upper end of smoothing function

        Returns:
        - step:         1d array of the step function
        '''
        step = np.zeros_like(lam)
        in_range = (lam >= lam_left) & (lam <= lam_right)
        x = (lam[in_range] - lam_left) / (lam_right - lam_left)
        step[in_range] = 6 * x**5 - 15 * x**4 + 10 * x**3
        step[lam > lam_right] = 1.0
        return step

    @classmethod
    def exp_bathymetry(cls, lam, lam_left, lam_right, width, slope, dist_taper):
        '''
        Compute tapered exponential slopes to define the bathymetry.

        Parameters:
        - lam:          coordinate value (e.g. lon or lat) 
        - lam_left:     lower end of lam
        - lam_right:    upper end of lam
        - width:        width of the basin in degrees longitude
        - slope:        slope of the exponential in degrees longitude
        - dist_taper:   distance from boundary where Hmax should be reached

        Returns:
        -exp:           value of the exponential function at `lam`
        '''
        exp = np.zeros_like(lam)

        norm = 1 + np.exp(-width / slope)

        left_taper = (lam >= lam_left) & (lam <= lam_left + dist_taper)
        step_left = cls.smooth_step(lam[left_taper], lam_left, lam_left + dist_taper)
        exp[left_taper] = (1 - (np.exp(-(lam[left_taper] - lam_left) / slope)) / norm) * (1 - step_left) + step_left

        right_taper = (lam >= lam_right - dist_taper) & (lam <= lam_right)
        step_right = 1 - cls.smooth_step(lam[right_taper], lam_right - dist_taper, lam_right)
        exp[right_taper] = (1 - (np.exp((lam[right_taper] - lam_right) / slope)) / norm) * (1 - step_right) + step_right

        main_basin = (lam > lam_left + dist_taper) & (lam < lam_right - dist_taper)
        exp[main_basin] = 1.0

        return exp
    
    @staticmethod
    def scurve(x, x0, dx):
        '''
        Returns 0 for x<x0 or x>x+dx, and a cubic in between.
        
        Taken from
        https://github.com/ocean-eddy-cpt/NeverWorld2/blob/main/docs/Wind%20profile.ipynb
        used in the Neverworld2 configuration. (Marques et al. 2022)
        '''
        s = np.minimum(1, np.maximum(0, (x-x0)/dx))
        return (3 - 2*s)*( s*s )

    @classmethod
    def tau_from_nodes(cls, phi_node, tau_node, phi):
        '''
        Returns a profile tau(phi) that uses s-curves between nodes.
        
        Taken from
        https://github.com/ocean-eddy-cpt/NeverWorld2/blob/main/docs/Wind%20profile.ipynb
        used in the Neverworld2 configuration. (Marques et al. 2022)
        '''
        tau = 0.*phi
        ks = 0
        for i in range(len(phi)):
            phi_i = phi[i]
            if phi_i>=phi_node[ks+1]:
                ks=min(len(phi_node)-2,ks+1)
            tau[i] = tau_node[ks] + ( tau_node[ks+1] - tau_node[ks]) * cls.scurve(phi_i, phi_node[ks], phi_node[ks+1]-phi_node[ks])
        return tau
    
    @staticmethod
    def seasonal_cycle(day_of_year=0.0, lag=0):
        '''
        Compute seasonal cycle.

        Parameters:
        - day_of_year:  days after dezember 31st [days]
        - lag:          shift the seasonal cycle [days]

        returns: seasonal_cycle as cosine(day_of_year)
        '''
        time_max = 5. * 30. + 21. + lag     # 21th june     at 24h in hours
        time_min = 11. * 30. + 21 + lag     # 21th december        in hours

        # for broadcasting
        day_of_year = np.asarray(day_of_year)
        
        seasonal_cycle = np.cos( (day_of_year - time_max) / (time_min - time_max) * np.pi )
        return(seasonal_cycle)
    

    
