'''
Created on Nov 29, 2018

@author: bry
'''
import sys
from collections import OrderedDict
import scipy.interpolate
import numpy as np

import includes.misc

def compute_rotation_matrix(angles):
    '''
    Generates the rotation matrix based on the angles supplied.
    '''
    # Assign the angles and convert to radians
    # Also called alpha angle and roughly corresponding to pattern phi angle
    yaw = angles[0] * np.pi / 180
    # Also called beta angle and roughly corresponding to pattern theta angle
    pitch = angles[1] * np.pi / 180
    # Also called psi angle. Does not map to any pattern angle
    roll = angles[2] * np.pi / 180
    # Make rotation matrixes
    rot_x = np.matrix([[1, 0, 0],
            [0, np.cos(yaw), np.sin(yaw)],
            [0, -np.sin(yaw), np.cos(yaw)]])
    rot_y = np.matrix([[np.cos(pitch), 0, -np.sin(pitch)],
            [0, 1, 0],
            [np.sin(pitch), 0, np.cos(pitch)]])
    rot_z = np.matrix([[np.cos(roll), np.sin(roll), 0],
            [-np.sin(roll), np.cos(roll), 0],
            [0, 0, 1]])
    # Total rotation
    rot = rot_x * rot_y * rot_z
    return rot

def rotate(far_field, angles):
    '''
    Warning:
        This is a private method and is NOT intended for external use.

    This static rotation method uses geometric definition of alpha, beta, psi angles or yaw,
    pitch roll angles correspondingly. The fixed XYZ coordinate system is the same used to
    define the radiation pattern in spherical coordinates.

    Arguments:
        far_field (array): pattern in the default format as defined for the ``far_field`` \
            matrix.
        angles (int,int,int): an interger triplet describing the three Euler angles of rotation.
            The order is --> [yaw, pitch, roll] as defined in [IEEESTDTest]_.

    Returns:
        far_field (array): rotated ``far_field``
    '''
    if angles is None:
        angles = [0,0,0]

    if not np.sum(angles) == 0:

        rot = compute_rotation_matrix(angles)

        # Make the new Etheta and Ephi arrays and determine the shape
        no_theta_samp = far_field.shape[1]
        no_phi_samp = far_field.shape[2]
        theta_step = int(180//(no_theta_samp-1))
        phi_step = int(360//(no_phi_samp-1))
        v_th = np.array(range(0,180+theta_step,theta_step))*np.pi/180
        v_ph = np.array(range(0,360,phi_step))*np.pi/180
        ff_rot = np.zeros(far_field.shape, dtype=complex)

        # Ittterate over each pattern value
        for i_2 in range(0,no_theta_samp):
            for j_2 in range(0,no_phi_samp):
                # Theta and phi spherical coordinates in the new pattern
                theta_2 = v_th[i_2]
                phi_2 = v_ph[j_2]

                # Find spherical coordinate unit vectors
                r_2  = np.matrix([  [ np.sin(theta_2)*np.cos(phi_2) ] , \
                                    [ np.sin(theta_2)*np.sin(phi_2) ] , \
                                    [     np.cos(theta_2)           ] ])# Defines column vector
                th_2 = np.matrix([  [ np.cos(theta_2)*np.cos(phi_2) ] , \
                                    [ np.cos(theta_2)*np.sin(phi_2) ] , \
                                    [ -np.sin(theta_2)              ] ])# Defines column vector
                ph_2 = np.matrix([  [          -np.sin(phi_2)       ] , \
                                    [            np.cos(phi_2)      ] , \
                                    [             0                 ] ])# Defines column vector

                # rotate them back into the old pattern coordinate system
                r_1  = rot * r_2
                th_1 = rot * th_2
                ph_1 = rot * ph_2

                # find coordinates in the old pattern
                theta_1 = float(np.arccos( r_1[2] ))
                phi_1   = float(np.arctan2( r_1[1], r_1[0] ))
                if phi_1 < 0:
                    phi_1 = phi_1 + 2*np.pi # we want only positive angles

                # find nearest neighbor and obtain indices
                i_1 = int(round( theta_1 / (theta_step*np.pi/180) ))
                j_1 = int(round( phi_1   / (phi_step*np.pi/180)   ))
                if j_1 == no_phi_samp:
                    j_1 = 0 # wraparound

                # special case at "north" and "south" "poles": j_1 can be inaccurate
                if i_1 == 0 or i_1 == no_theta_samp-1:
                    j_1 = 0
                # theta and phi coordinates in the old pattern
                theta_1 = v_th[i_1]
                phi_1   = v_ph[j_1]

                # unit vectors in the old pattern defined as column vectors
                th1_orig = np.matrix([  [np.cos(theta_1)*np.cos(phi_1)] , \
                                        [np.cos(theta_1)*np.sin(phi_1)] , \
                                        [    -np.sin(theta_1)         ]  ])
                ph1_orig = np.matrix([  [          -np.sin(phi_1)     ] , \
                                        [           np.cos(phi_1)     ] , \
                                        [           0                 ]  ])

                # angle cosines
                cos_th_th = float(th_1.T * th1_orig)
                cos_th_ph = float(th_1.T * ph1_orig)
                cos_ph_th = float(ph_1.T * th1_orig)
                cos_ph_ph = float(ph_1.T * ph1_orig)

                #  fill the new Etheta2 and Ephi2 arrays
                ff_rot[0,i_2,j_2] = \
                    far_field[0,i_1,j_1] * cos_th_th + far_field[1,i_1,j_1] * cos_th_ph # Theta
                ff_rot[1,i_2,j_2] = \
                    far_field[0,i_1,j_1] * cos_ph_th + far_field[1,i_1,j_1] * cos_ph_ph #Phi
    else:
        ff_rot = far_field

    return ff_rot

def convert_field_format(field, field_format, temp_field):
    '''
    Warning:
        This is a private method and is NOT intended for external use.

    Arguments:
        field (str): the type of field requested. See ``__convert_field`` for list of \
            possible values.

        field_format (str): The format of the output field. Following definitions are according
            to the Touchstone file format:

            - **dB** - for logarithmic representation of the magnitude and degrees angle for \
                the phase (where applicable). This is the default format.

            - **MA** - for linear representation of the magnitude and degrees angle for the\
                phase (where applicable)

            - **RI** - for complex real and imaginary representaion of the field \
                (where applicable)

    Returns:
        list: list containing a two dimensional array with the theta/phi values of the \
            field along the columns and rows correspondingly. The number columns/rows \
            corresponds to the pattern resolution along the theta and phi angles.
    '''
    if field_format == 'dB':
        temp_field = includes.misc.lin2log(field, temp_field)
    elif field_format == 'MA' and any(field in s for s in ['AR', 'Eth', 'Eph', 'EL', 'ER']):
        temp_field = np.abs(temp_field)
    elif field_format == 'RI' and any(field in s for s in ['Gabs', 'Gth', 'Gph', 'GL', 'GR', \
                                                           'Eabs', 'AR']):
        sys.exit('RI format is ONLY supported for individual electric fields')
    return temp_field

def convert_field(field, field_format, pattern):
    '''
    .. _func_convert_field:

    Converst the input ``pattern`` to one of the specified ``field`` using the ``field_format``.

    Arguments:
        field (str): The field format requested. Options are:

            - **Gabs** - Total Realized Gain is the sum of the two polarizations. It can be in \
            linear or logarithmic magnitude. Phase is NOT defined. Real/Imaginary output is \
            NOT defined. This is the default output.

            - **Gth** - Theta Realized Gain is for linear polarization in spherical coordinate \
            system. It can be in linear or logarithmic magnitude. The phase is NOT defined \
            for the gain but rather for the associated E field. The phase is returned in \
            degrees. Real/Imaginary output is NOT defined.

            - **Gph** - Phi Realized Gain is for linear polarization in spherical coordinate \
            system. It can be in linear or logarithmic magnitude. The phase is NOT defined \
            for the gain but rather for the associated E field. The phase is returned in \
            degrees. Real/Imaginary output is NOT defined.

            - **GL** - Left Realized Gain is for circular polarization in spherical coordinate \
            system. It can be in linear or logarithmic magnitude. The phase is NOT defined \
            for the gain but rather for the associated E field. The phase is returned in \
            degrees. Real/Imaginary output is NOT defined.

            - **GR** - Right Realized Gain is for circular polarization in spherical \
            coordinate system. It can be in linear or logarithmic magnitude. The phase is NOT \
            defined for the gain but rather for the associated E field. The phase is returned \
            in degrees. Real/Imaginary output is NOT defined.

            - **AR** - Axial Ratio is mostly for circularly polarized antennas. It can be in \
            linear or logarithmic magnitude. Phase is NOT defined. Real/Imaginary output is \
            NOT defined.

            - **Eabs** - total electric field is the sum of the two polarizations. It can be \
            in linear or logarithmic magnitude. When in logarithmic scale this is the same as \
            the corresponding Gabs. Phase is NOT defined. Real/Imaginary output is NOT \
            defined.

            - **Eth** - theta electric field is for linear polarization in spherical coordinate\
            system. It can be as a complex number (real/imaginary output) or as linear or \
            logarithmic magnitude with phase. The phase is returned in degrees.

            - **Eph** - phi electric field is for linear polarization in spherical coordinate \
            system. It can be as a complex number (real/imaginary output) or as linear or \
            logarithmic magnitude with phase. The phase is returned in degrees.

            - **EL** - left electric field is for circular polarization in spherical coordinate\
            system. It can be as a complex number (real/imaginary output) or as linear or \
            logarithmic magnitude with phase. The phase is returned in degrees.

            - **ER** - right electric field is for circular polarization in spherical \
            coordinate system. It can be as a complex number (real/imaginary output) or as \
            linear or logarithmic magnitude with phase. The phase is returned in degrees.

        pattern (int,int,int): the patttern from the local database to be converted. See \
            ``__get_raw_field_data`` for details on expected input.

    Returns:
        tuple: tuple containing:

            field (list): a two dimensional array containing the theta/phi values of the field \
                along the columns and rows correspondingly. The number columns/rows corresponds\
                to the pattern resolution along the theta and phi angles.

            phase (list): a two dimensional array containing the theta/phi values of the field\
                phase along the columns and rows correspondingly. The number columns/rows\
                corresponds to the pattern resolution along the theta and phi angles. See field\
                definitions above for when phase is defined. When phase is NOT defined an empty\
                list is returned.
    '''

    e_left = (1 / np.sqrt(2)) * (pattern[0] - 1j * pattern[1]) # converts to left E field
    e_right = (1 / np.sqrt(2)) * (pattern[0] + 1j * pattern[1]) # converts to right E field
    if field == 'Gabs':
        temp_field = np.abs(pattern[0]) ** 2 + np.abs(pattern[1]) ** 2
        temp_phase = [] # phase is NOT defined for this field
    elif field == 'Gth':
        temp_field = np.abs(pattern[0]) ** 2
        temp_phase = np.angle(pattern[0], deg=1)
    elif field == 'Gph':
        temp_field = np.abs(pattern[1]) ** 2
        temp_phase = np.angle(pattern[1], deg=1)
    elif field == 'GL':
        temp_field = np.abs(e_left) ** 2
        temp_phase = np.angle(e_left, deg=1)
    elif field == 'GR':
        temp_field = np.abs(e_right) ** 2
        temp_phase = np.angle(e_right, deg=1)
    elif field == 'AR':
        temp_field = (np.abs(e_right) + np.abs(e_left)) / (np.abs(e_right) - np.abs(e_left))
        temp_phase = [] # phase is NOT defined for this field
    elif field == 'Eabs':
        temp_field = np.abs(pattern[0]) ** 2 + np.abs(pattern[1]) ** 2
        temp_phase = [] # phase is NOT defined for this field
    elif field == 'Eth':
        temp_field = pattern[0]
        temp_phase = np.angle(pattern[0], deg=1)
    elif field == 'Eph':
        temp_field = pattern[1]
        temp_phase = np.angle(pattern[1], deg=1)
    elif field == 'EL':
        temp_field = e_left
        temp_phase = np.angle(e_left, deg=1)
    elif field == 'ER':
        temp_field = e_right
        temp_phase = np.angle(e_right, deg=1)

    temp_field = includes.misc.fix_singularities(temp_field)
    temp_field = convert_field_format(field, field_format, temp_field)

    return temp_field, temp_phase

def get_far_field_components(far_field, method = 'phasor'):
    '''
    Returns Magnitude/Angle or Real/Imaginary component matrixes for both Theta/Phi angles.
    '''
    components = []
    if method == 'phasor':
        # Extract phasors - magnitude and phase
        components.extend(convert_field('Eth', 'MA', far_field))
        components.extend(convert_field('Eph', 'MA', far_field))
        components[1] = components[1] * np.pi / 180
        components[3] = components[3] * np.pi / 180
    elif method == 'RI':
        # in RI format
        components.append(np.real(far_field[0]))
        components.append(np.imag(far_field[0]))
        components.append(np.real(far_field[1]))
        components.append(np.imag(far_field[1]))
    else:
        sys.exit('Unknown complex interpolation method! --> {}'.format(method))
    return components

def combine_far_field_components(components, method = 'phasor'):
    '''
    Reconstructs the electric fields in the patterns.py format. This operation is the reverse of
    ``get_far_field_components``
    '''
    if method == 'phasor':
        far_field = np.stack((
                                components[0] * np.cos(components[1]) + \
                                (1j * components[0] * np.sin(components[1])),\
                                components[2] * np.cos(components[3]) + \
                                (1j * components[2] * np.sin(components[3]))))
    elif method == 'RI':
        far_field = np.stack((
            components[0] + (1j * components[1]),
            components[2] + (1j * components[3])))
    return far_field

def corelation_calculation(ff_1, ff_2, env, xpd): # ff_1/ff_2/env shapes MUST be the same!!!!
    '''
    Computes complex radiation pattern correlation based on J.B Andersen's paper from '87.

    Todo:
        add more detail on inputs and output and provide IEEE reference.
    '''
    resolution = includes.misc.get_angle_steps(ff_1[0], sphere_type = 'closed')
    vth = np.array(range(0,int(180+resolution[0]),int(resolution[0])))*np.pi/180
    # vph is defined 0:360 for the integration only. For normal matrix storage the 360 degrees
    # is not needed
    vph = np.array(range(0, int(360+resolution[1]), int(resolution[1])))*np.pi/180

    numerator = np.zeros(ff_1[0].shape, dtype=complex)
    denominator1 = np.zeros(ff_1[0].shape, dtype=complex)
    denominator2 = np.zeros(ff_1[0].shape, dtype=complex)

    for i in range(vth.__len__()):
        numerator[i,:] = (xpd * ff_1[0][i,:]*np.conj( ff_2[0][i,:])*env[0][i,:] +     \
                                ff_1[1][i,:]*np.conj( ff_2[1][i,:])*env[1][i,:]   ) * \
                                np.sin(vth[i])
        denominator1[i,:] = (xpd *  ff_1[0][i,:]*np.conj( ff_1[0][i,:])*env[0][i,:]  +    \
                                    ff_1[1][i,:]*np.conj( ff_1[1][i,:])*env[1][i,:]   ) * \
                                    np.sin(vth[i])
        denominator2[i,:] = (xpd *  ff_2[0][i,:]*np.conj( ff_2[0][i,:])*env[0][i,:]  +    \
                                    ff_2[1][i,:]*np.conj( ff_2[1][i,:])*env[1][i,:]   ) * \
                                    np.sin(vth[i])

    # integrate
    numerator = np.trapz(np.trapz(numerator, vth, axis=0), vph)
    denominator1 = np.trapz(np.trapz(denominator1, vth, axis=0), vph)
    denominator2 = np.trapz(np.trapz(denominator2, vth, axis=0), vph)

    corr = (numerator)/np.sqrt(denominator1*denominator2)

    return corr

def antenna_temperature_calculation(gain, brightness_temperature):
    '''
    Computes the antenna temperature :math:`T_A` as defined in [BALANIS]_.

    .. math::
        T_A = \\frac{\\int_{0}^{2\\pi} \\int_{0}^{\\pi} T_B( \\theta, \\phi)G( \\theta,\\phi) \
        sin( \\theta) d\\theta d\\phi}{\\int_{0}^{2\\pi} \\int_{0}^{\\pi} \
        G( \\theta,\\phi)sin( \\theta) d\\theta d\\phi}

    Args:
        gain (array): numpy array holding the ``Gabs`` radiation pattern (:math:`G(\\theta,\\phi)`)
            of the antenna - see :func:`~includes.field_operations.convert_field`.
        brightness_temperature (array): temperature pattern (:math:`T_B`) generated with
            :func:`~includes.field_operations.generate_brightness_temperature`.

    Returns:
        float: the calculated antenna temperature. **NOTE** the returned temeperature does NOT
        include the contributions of the physical antenna temperature.
    '''
    vth,vph = includes.misc.get_angle_vectors(gain, sphere_type = 'open', return_mesh = False)

    numerator = np.zeros(gain.shape, dtype=float)
    denominator = np.zeros(gain.shape, dtype=float)

    for i in range(vth.__len__()):
        numerator[i,:] = brightness_temperature[i,:] * gain[i,:] * np.sin(vth[i])
        denominator[i,:] = gain[i,:] * np.sin(vth[i])
    # integrate
    numerator = np.trapz(np.trapz(numerator, vth, axis=0), vph)
    denominator = np.trapz(np.trapz(denominator, vth, axis=0), vph)

    antenna_temperature = numerator / denominator

    return antenna_temperature

def generate_brightness_temperature(model = None, shape = None):
    '''
    This method generates a 3D brightness temperature specified by the ``model`` and dimenssioned
    according to ``shape``. The ``shape`` should be such that the generated brightness temperature
    can be directly used to multiple an antenna pattern and calculate antenna noise temperature with
    :func:`~includes.field_operations.antenna_temperature_calculation`.

    Args:
        model (dict): a python dictionary describing the brightness temperature model. The options\
            are:

                - type (:obj:`str`): a string describing the model type. Currently this can be\
                    *'uniform'* (*default*) or 'cone'.

                    - *'uniform'* brightness model generates same ``temperature`` from all \
                        directions.

                    - 'cone' brightness model generates a cone with values ``temperature`` in a \
                        ``diameter`` around ``direction``. Values outside the cone are set to \
                        ``cosmic background``

                - temperature (:obj:`int`): the value of the desired brightness temperature. By \
                    default this is set equal to the ``cosmic background`` of 5 dgrees Kelvin.

                - cosmic background (:obj:`int`): the value of the cosmic background temperature.\
                    By default this is set to 5 dgrees Kelvin [NEEDCITATION]_.

                - diameter (:obj:`int`): the diameter of the ``cone`` type model expressed in \
                    degrees around ``direction``.

                - direction (:obj:`int`, :obj:`int`): direction of the cone on a sphere defined as \
                    Theta and Phi angles of the main direction.
        shape (list): a ist of two values (theta,phi) for the size of the array to be generated - \
            this typically comes from the pattern of the antenna.
    Returns:
        array: a numpy array of ``shape`` holding the brightness temperature for each theta/phi
        andgle
    '''
    # State the defaults
    model_defaults = {  'type':'uniform',
                    'temperature':5,
                    'cosmic background':5,
                    'diameter':20,
                    'direction':[0,0]
                    }
    model = includes.misc.update_params(model_defaults,model)
    if shape is None:
        shape = [37,73] # corresponds to 5 degrees resolution for both theta and phi

    if model['type'] == 'uniform':
        t_brightness =  np.tile(model['temperature'], shape)
    elif model['type'] == 'cone':
        resolution = includes.misc.get_angle_steps_from_shape(shape, sphere_type='open')
        vth = np.tile(model['cosmic background'], shape[0])
        vth[0:int(model['diameter']/resolution[0])] = model['temperature']
        t_brightness = np.tile(vth, shape[1])
        t_brightness = np.reshape(t_brightness, shape ,'F')
        includes.misc.check_angles_input(model['direction'], 'theta/phi direction')
        t_brightness = np.roll(t_brightness, model['direction'],axis=(0,1))

    return t_brightness

def generate_environment(environment = None, shape = None):
    '''
    Generates a power distribution function for an environemnt. Currently only isotropic
    environment is supported with values of 1/4*pi for each theta and phi polarizations.
    '''
    if environment is None:
        env = 'isotropic'
        env_params = OrderedDict()
        env_params['XPD'] = 1
    else:
        env = environment[0]
        env_params = environment[1]
    if shape is None:
        shape = [37,73] # corresponds to 5 degrees resolution for both theta and phi

    if env == 'isotropic':
        p_th = np.tile(1/(4*np.pi), shape)
        env_pdp = np.stack((p_th,p_th))

    return env_pdp, env_params['XPD']

def interpolate_field(interp_v, matrix,pat_resolution):
    '''
    Interpolates a ``matrix`` of certain ``resolution`` to the sizes in ``interp_v``. This function
    is the lowest level and works on a single matrix of real numbers.
    '''
    matrix = includes.misc.close_sphere(matrix)

    # Extract the resolution and generate stepping vectors for a closed sphere
    vth = np.array(range(0, 180+pat_resolution[0], pat_resolution[0]))*np.pi/180
    vph = np.array(range(0, 360+pat_resolution[1], pat_resolution[1]))*np.pi/180

    no_patt_samples = interp_v[0].shape

    lut_matrix = scipy.interpolate.RectBivariateSpline(vth,vph,matrix,kx=1,ky=1,s=0)
    matrix = lut_matrix.ev(interp_v[0].ravel(),interp_v[1].ravel())\
        .reshape((no_patt_samples[0], no_patt_samples[1])).T

    return matrix

def interpolate_complex_field(interp_v, far_field,pat_resolution, method = 'phasor'):
    '''
    Special function handling interpolation of complex electric fields matrixes by decomposing to
    individual real matrix components and re-assembling after the interpolation.
    '''
    components = get_far_field_components(far_field, method)

    for idx, f_to_interp in enumerate(components):
        components[idx] = interpolate_field(interp_v, f_to_interp, pat_resolution)

    far_field = combine_far_field_components(components, method)

    return far_field
