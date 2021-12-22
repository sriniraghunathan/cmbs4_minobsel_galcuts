def get_lat_based_masks(nside, max_lat_mask = 50., delta_lat = 10., rotate_to_radec = True):
    """
    #get lat-based masks in ra/dec coordinates
    """
    lat_arr = np.arange( 0., max_lat_mask+1., delta_lat)

    hmask_dic = {}
    npix = H.nside2npix( nside )
    phi_deg, theta_deg = H.pix2ang( nside, np.arange(npix), lonlat = 1 )
    
    for lval in lat_arr:
        curr_mask = np.ones( npix )
        unmask_pixels = np.where( (theta_deg>=-lval) & (theta_deg<=lval) )[0]
        curr_mask[unmask_pixels] = 0.

        if rotate_to_radec:#rotate to ra/dec
            curr_mask = healpix_rotate_coords(curr_mask, coord = ['G', 'C'])

        '''
        H.mollview( curr_mask, coord = ['C', 'G'], sub = (2,2,1) ); H.graticule();
        H.mollview( curr_mask, sub = (2,2,2) ); H.graticule(); show(); #sys.exit()
        '''
        hmask_dic[lval] = curr_mask

    return hmask_dic #these are lat masks in ra/dec coordinates

def get_lat_based_masks_v1(nside, min_lat = -90., max_lat = 45., delta_lat = 15., rotate_to_radec = True, use_lat_steps = False):
    """
    #get lat-based masks in ra/dec coordinates
    """
    lat_arr = np.arange( min_lat, max_lat + 1., delta_lat )

    hmask_dic = {}
    npix = H.nside2npix( nside )
    phi_deg, theta_deg = H.pix2ang( nside, np.arange(npix), lonlat = 1 )
    
    lat_mask_arr = []
    for l1 in lat_arr[:-1]:
        l2 = l1 + delta_lat
        if not use_lat_steps:
            l1 = lat_arr[0] #force l1 to be -90 deg.
        curr_mask = np.zeros( npix )
        unmask_pixels = np.where( (theta_deg>=l1) & (theta_deg<l2) )[0]
        curr_mask[unmask_pixels] = 1.

        if rotate_to_radec:#rotate to ra/dec
            curr_mask = healpix_rotate_coords(curr_mask, coord = ['G', 'C'])

        #H.mollview( curr_mask, sub = (2,2,1)); H.graticule();
        #H.mollview( cmbs4_hit_map, sub = (2,2,2)); H.graticule(); show()

        #H.mollview( curr_mask ); H.graticule(); show()
        print((l1, l2))
        hmask_dic[(l1, l2)] = curr_mask

        lat_mask_arr.append(curr_mask)

    if use_lat_steps: #finally add all masks together
        curr_mask = np.sum(lat_mask_arr, axis = 0)
        curr_mask[curr_mask>threshold] = 1.
        curr_mask[curr_mask<threshold] = 0.
        lat_mask_arr.append( curr_mask )
        l1, l2 = lat_arr[0], lat_arr[-1]
        hmask_dic[(l1, l2)] = curr_mask
        
    #H.mollview( curr_mask ); H.graticule(); show(); sys.exit()

    #return np.asarray( lat_mask_arr ) #these are lat masks in ra/dec coordinates
    return hmask_dic #these are lat masks in ra/dec coordinates
def healpix_rotate_coords(hmap, coord):
    """
    coord = ['C', 'G'] to convert a map in radec to galactic coordinates or vice versa.
    """

    #get map pixel
    pixel = np.arange(len(hmap))

    #get angles in this map first
    nside = H.get_nside(hmap)
    angles = H.pix2ang(nside, pixel)

    #roate the angles to the desired new coordinate
    rotated_angles = H.Rotator(coord=coord)(*angles)

    #get the rotated pixel values
    rotated_pixel = H.ang2pix(nside, *rotated_angles)

    #initialise new map
    rot_hmap = np.zeros(len(pixel))

    #push the original map pixel to the new map (in the rotated pixel positions)
    rot_hmap[rotated_pixel] = hmap[pixel]
    if (1):
        rot_hmap = H.smoothing(rot_hmap, fwhm = np.radians(beam_to_use_for_smoothing))
        rot_hmap[rot_hmap>threshold] = 1.
        rot_hmap[rot_hmap<threshold] = 0.

    return rot_hmap

def set_telescope_observer(lat = -89.991066, lon = -44.65, elevation = 2835.0, astropy = 1):#, location = 'southpole'):
    '''
    if location == 'southpole':
        lat, lon, elevation = -22.57, -67.47, 5200. #Chile
    elif location == 'chile':
        lat, lon, elevation = -22.57, -67.47, 5200. #Chile
    '''
    if not astropy:
        splat = ephem.Observer()
        splat.lat = np.radians(lat)
        splat.lon = np.radians(lon)
        splat.elevation = elevation # in meters
    else:
        ##splat = EarthLocation.from_geodetic(lon* u.degree, lat* u.degree, height = elevation*u.meter)
        splat = coord.EarthLocation(lat=lat*u.deg,lon=lon*u.deg, height=elevation*u.meter)

    return splat

def get_hmask(nside, min_obs_el, max_obs_el = 85., timearr = None, ra_dec_arr = None, yyyy = 2029, delta_mm = 1, delta_dd = 30., delta_hh = 24., location = 'southpole'):#, epoch = 'J2000'):
    if ra_dec_arr is None:
        raarr = np.arange(-180, 180., 0.1)
        decarr = np.arange(-90., 90., 0.1)

        ra_dec_arr = []
        for r in raarr:
            for d in decarr:
                ra_dec_arr.append( [r, d])
        ra_dec_arr = np.asarray( ra_dec_arr )

    if timearr is None: #does not matter for pole.
        timearr = []
        for mm in np.arange(1, 13, delta_mm):
            for dd in np.arange(1, 30, delta_dd):
                for hh in np.arange(1, 24, delta_hh):
                    tval = '%04d-%02d-%02dT%02d:00:00' %(yyyy, mm, dd, hh)
                    try:
                        t = Time(tval, format='isot', scale='utc')
                        timearr.append( tval )
                    except:
                        pass
        timearr = Time(timearr, format='isot', scale='utc')

    astropy_coord_radec = SkyCoord(ra_dec_arr[:,0], ra_dec_arr[:,1], unit='deg', frame='icrs')
    splat = set_telescope_observer()
    npix = H.nside2npix(nside)
    hmask = np.zeros(npix)

    for t in timearr:
        astropy_coord_altaz = astropy_coord_radec.transform_to(coord.AltAz(obstime=t,location=splat))
        az, el = astropy_coord_altaz.az.value, astropy_coord_altaz.alt.value
        good_inds = np.where( (el>=min_obs_el) & (el<=max_obs_el) )[0]
        good_ra_dec = ra_dec_arr[good_inds]

        good_pixels = H.ang2pix(nside, np.radians(90.-good_ra_dec[:,1]), np.radians(good_ra_dec[:,0]))

        hmask[good_pixels] = 1.
        if location == 'southpole': break

    hmask = H.smoothing(hmask, fwhm = np.radians(beam_to_use_for_smoothing))    
    hmask[hmask<threshold] = 0.
    hmask[hmask>threshold] = 1.

    return hmask

###########################################################################
###########################################################################

from pylab import *
from astropy.time import Time
import astropy, numpy as np, os, sys
import healpy as H
from astropy import units as u
from astropy import coordinates as coord
from astropy.coordinates import SkyCoord, EarthLocation

#use_lat_steps = False #True #False
beam_to_use_for_smoothing = 0.5 #some rough smoothing
threshold = 0.001
nside = 512
lmax = 512
dust_145ghz_fname = 'data/cmbs4_dust_uKCMB_LAT-MFL2_nside512_dust_0000.fits'
dust_145ghz = H.read_map(dust_145ghz_fname, field = (0, 1, 2))
#dust_145ghz = H.ud_grade(dust_145ghz, nside)

lat_mask_dic = get_lat_based_masks(nside)#, use_lat_steps = use_lat_steps)
min_obs_el_err = [20., 25., 30., 35., 40.]
op_dic = {}
op_dic['hit_masks'] = {}
op_dic['lat_masks'] = lat_mask_dic
op_dic['cl'] = {}
for min_obs_el in min_obs_el_err:
    op_dic['cl'][min_obs_el] = {}
    hmask = get_hmask(nside, min_obs_el = min_obs_el)
    op_dic['hit_masks'][min_obs_el] = hmask
    fsky1 = len( np.where(hmask>0.)[0] )/len(hmask) #all inds with more than xx per cent hit    
    fsky = np.sum( hmask[hmask>0.] )/len(hmask) #all inds with more than xx per cent hit    
    print(beam_to_use_for_smoothing, 0.5*(1-np.sin(np.radians(min_obs_el))), min_obs_el, fsky1, fsky, np.mean(hmask)); #sys.exit()

    '''
    H.mollview(hmask, min = 0., max = 1., coord = ['C', 'G'], sub=(2,3,1)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); 
    H.mollview(lat_mask_dic[(-90.0, 45.0)], min = 0., max = 1., coord = ['C', 'G'], sub=(2,3,2)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); 
    H.mollview(hmask * lat_mask_dic[(-90.0, 45.0)], min = 0., max = 1., coord = ['C', 'G'], sub=(2,3,3)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky));

    H.mollview(hmask, min = 0., max = 1., sub=(2,3,4)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); 
    H.mollview(lat_mask_dic[(-90.0, 45.0)], min = 0., max = 1., sub=(2,3,5)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); 
    H.mollview(hmask * lat_mask_dic[(-90.0, 45.0)], min = 0., max = 1., sub=(2,3,6)); H.graticule(); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); show(); sys.exit()
    '''

    #for (l1, l2) in lat_mask_dic:
        #lat_mask = lat_mask_dic[(l1, l2)]
    for bval in lat_mask_dic:
        lat_mask = lat_mask_dic[bval]
        '''
        fsky = np.sum( np.where(hmask>0.)[0] )/len(hmask) #all inds with more than xx per cent hit    
        H.mollview(hmask, min = 0., max = 1., sub = (1,3,1)); H.graticule(verbose = False); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky));

        H.mollview(hmask_with_gal_cut, min = 0., max = 1., sub = (1,3,2)); H.graticule(verbose = False); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); 

        curr_dust_map = hmask * lat_mask * dust_145ghz
        H.mollview(curr_dust_map, sub = (1,3,3)); H.graticule(verbose = False); title(r'Min. el = %g; fsky = %.3f' %(min_obs_el, fsky)); show(); 
        H.mollview(curr_dust_map); H.graticule(verbose = False); show()
        '''
        hmask_with_gal_cut = hmask * lat_mask
        hmask_with_gal_cut[hmask_with_gal_cut<threshold] = 0.
        #fsky = len( np.where(hmask_with_gal_cut>0.)[0] )/len(hmask_with_gal_cut) #all inds with more than xx per cent hit
        fsky = np.sum( hmask_with_gal_cut[hmask_with_gal_cut>0.])/len(hmask_with_gal_cut) #all inds with more than xx per cent hit

        #print('\t(l1, l2) = (%g, %g) and fsky = %.3f' %(l1, l2, fsky))
        print('\tLast mask = +/- %g deg. and fsky = %.3f' %(bval, fsky))
        
        curr_mask = hmask * lat_mask 
        curr_dust_map = curr_mask * np.copy( dust_145ghz )
        cl = H.anafast(curr_dust_map, lmax = lmax)
        cl = cl/fsky
        #op_dic['cl'][min_obs_el][(l1,l2)] = cl
        op_dic['cl'][min_obs_el][bval] = cl

op_dic['dust_145ghz'] = dust_145ghz
'''
if not use_lat_steps:
    opfname = 'results/splat_hitmask_galmask_galdustcl.npy'
else:
    opfname = 'results/splat_hitmask_galmask_galdustcl_latsteps.npy'
'''
opfname = 'results/splat_hitmask_galmask_galdustcl.npy'
np.save(opfname, op_dic)
sys.exit()



