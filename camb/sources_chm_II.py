from .baseconfig import F2003Class, fortran_class, numpy_1d, np, numpy_1d_or_null
from ctypes import POINTER, c_int, c_double, byref
#from ctypes import c_int_or_null, c_double_or_null # CHM


class SourceWindow(F2003Class):
    """
    Abstract base class for a number count/lensing/21cm source window function.
    A list of instances of these classes can be assigned to the SourceWindows field of :class:`.model.CAMBparams`.

    Note that source windows can currently only be used in flat models.
    """
    _fields_ = [("source_type", c_int, {"names": ["21cm", "counts", "lensing", "arf", "vlos"], "start": 1}),
                ("bias", c_double),
                ("dlog10Ndm", c_double)]

    _fortran_class_module_ = "SourceWindows"
    _fortran_class_name_ = "TSourceWindow"


@fortran_class
class GaussianSourceWindow(SourceWindow):
    """
    A Gaussian W(z) source window function.
    """
    _fields_ = [("redshift", c_double),
                ("sigma", c_double)]

    _fortran_class_name_ = "TGaussianSourceWindow"


@fortran_class
class SplinedSourceWindow(SourceWindow):
    """
    A numerical W(z) source window function constructed by interpolation from a numerical table.
    """
    _fortran_class_name_ = "TSplinedSourceWindow"

    #_methods_ = [("SetTable", [POINTER(c_int), numpy_1d, numpy_1d, numpy_1d_or_null ])] 
    _methods_ = [("SetTable", [POINTER(c_int), numpy_1d, numpy_1d, numpy_1d, numpy_1d, numpy_1d_or_null ])] # chm

    def __init__(self, **kwargs):
        z = kwargs.pop('z', None)
        if z is not None:
            #self.set_table(z, kwargs.pop('W'), kwargs.pop('bias_z', None))
            self.set_table(z, kwargs.pop('W'), kwargs.pop('bias_z', None), kwargs.pop('sigma_Errz', None), \
                kwargs.pop('i_Lorentz', None)) # chm, added "sigma_Errz" and "i_Lorentz"
        super().__init__(**kwargs)

    def set_table(self, z, W, bias_z=None, sigma_Errz=None, i_Lorentz=None): # chm, added "sigma_Errz" and "i_Lorentz"
    #def set_table(self, z, W, bias_z=None): 
        """
        Set arrays of z and W(z) for cublic spline interpolation. Note that W(z) is the total count distribution
        observed, not a fractional selection function on an underlying distribution.

        :param z: array of redshift values (monotonically increasing)
        :param W: array of window function values. It must be well enough sampled to smoothly cubic-spline interpolate
        :param bias_z: optional array of bias values at each z
        - CHM: 
        :param sigma_Errz: optional double scalar introduce photo-z error RMS in units of (1+z), i.e, 
                                                                            sigma_Errz*(1+z)=sigma_Errz_measured
        :param i_Lorentz: optional integer, if ==1 activates Lorentzian photo-z error PDF, otherwise a Gaussian in assumed
        - CHM ____
        """
        if len(W) != len(z) or z[-1] < z[1] or len(z) < 5:
            raise ValueError(
                "Redshifts must be well sampled and in ascending order, with window function the same length as z")
        if bias_z is not None:
            bias_z = np.asarray(bias_z, dtype=np.float64)
            if len(bias_z) != len(z):
                raise ValueError("bias array must be same size as the redshift array")
        # CHM
        if sigma_Errz is not None:
            sigma_Errz = np.asarray(sigma_Errz,dtype=np.float64)
        else:
            sigma_Errz=np.asarray(1e-5,dtype=np.float64)
        if i_Lorentz is not None:
            i_Lorentz = np.asarray(i_Lorentz,dtype=np.float64)
        else:
            i_Lorentz =  np.asarray(0,dtype=np.float64)
        #
        #self.f_SetTable(byref(c_int(len(z))), np.asarray(z, dtype=np.float64), np.asarray(W, dtype=np.float64), bias_z)
        self.f_SetTable(byref(c_int(len(z))), np.asarray(z, dtype=np.float64), np.asarray(W, dtype=np.float64),  \
                    sigma_Errz, i_Lorentz, bias_z )
        # CHM ____
