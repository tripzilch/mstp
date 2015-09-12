from pylab import *
from scipy.misc import imsave, toimage
from scipy import ndimage as nd
from scipy.interpolate import interp1d
from StringIO import StringIO
from functools import partial
import sys
import re

def box3(x, sigma, output, mode='wrap'):
    """Approximation of Gaussian filter by 3 box blurs.

    x -- input array
    sigma -- approx standard deviation

    The relationship between sigma and "n" (the number used to
    determine the box blur sizes) is as follows:
    
    sigma(n) = 0.345586 * n + 0.525521
    n(sigma) = 2.893638 * sigma - 1.520667

    (as empirically determined by finding the sigma for which
    scipy.ndimage.gaussian_filter gave the smallest max squared
    error.)
    """
    n = int(round(2.893638 * sigma - 1.520667))
    a, b, c = 1 + 2 * (n / 3), 1 + 2 * ((n + 1) / 3), 1 + 2 * ((n + 2) / 3)
    nd.uniform_filter(x, a, mode=mode, output=output) 
    nd.uniform_filter(output, b, mode=mode, output=output)
    nd.uniform_filter(output, c, mode=mode, output=output)
    
def boxG(x, sigma, output, mode='wrap'):
    """Gaussian filter if sigma < 5 else box3 filter."""
    f = nd.gaussian_filter if sigma < 5 else box3
    return f(x, sigma, output=output, mode=mode)

def box2(x, sigma, output, mode='wrap'):
    """Approximation of Gaussian filter by 2 box blurs.
    
    sigma(n) = 0.436589 * n + 0.441425
    n(sigma) = 2.290485 * sigma - 1.011077
    
    See box3 documentation for further details.    
    """
    n = int(round(2.290485 * sigma - 1.011077))
    a, b = 1 + 2 * (n / 2), 1 + 2 * ((n + 1) / 2)
    nd.uniform_filter(x, a, mode=mode, output=output) 
    nd.uniform_filter(output, b, mode=mode, output=output)
    
def symmetry_xform(z, N=3, mode='wrap'):
    """Blends input with random rotated version for N-fold symmetry."""
    return (z + nd.rotate(z, randint(1, N) * 360 / N, reshape=False, mode=mode, order=1)) * 0.5
    
    
    
class MSTP(object):
    """Class for representing and running Multiscale Turing Patterns.

    See [lnk] and [link] for explanation.
    """
    
    def __init__(self,
            shape=[300,300],
            ra=[ 64, 24,  9,  3,  1],
            ri=[ 96, 32, 12,4.5,1.5],
            dt=[.04,.03,.03,.02,.02], # wt=[1],
            pal=[[1,1,0],[1,0,0],[1,0,1],[0,0,1],[0,1,1],[1,1,1]]):
        """Create a new multiscale turing pattern.

        Keyword arguments:
        ra -- list of activator radii for each scale
        ri -- list of inhibitor radii for each scale
        dt -- list of timestep per for each scale
#        wt -- weight per scale (default None is equal weights)
        shape -- shape of buffer, initialized with uniform noise
        """

        "Greyscale buffer that contains the actual Multiscale Turing Pattern."
        self.z = rand(*shape)
        "Colour buffer with RGB values tracking the colour of local scales."
        self.c = ones(list(shape)+[3])
        "Timestep per scale."
        self.dt = array(dt)
#        "Weight per scale."
#        self.wt = array(wt)
        "Activator radii"
        self.ra = ra
        "Inhibitor radii"
        self.ri = ri
        "Colourmap of scale to RGB."
        self.pal = array(pal)
        "Transform function before filter."
        self._xform = lambda z: z
        "Filter function."
        self._filter = boxG
        "Colour buffer update speed."
        self._dc = .04 
        # init these as instance variables so they don't have to be
        # allocated on each call to self.step()
        self._tmp = zeros_like(self.z)
        self._min_var = zeros_like(self.z).astype(int)
        self._variance = zeros([len(ra)] + list(shape))
        
    def reset(self):
        """Re-initializes greyscale and colour buffers."""
        self.z = rand(*self.z.shape)
        self.c = ones_like(self.c)

    def step(self, speed=1.0):
        """Advance the multiscale Turing Pattern by one step."""
        # shortcut for slice indexing
        _XY = s_[:,newaxis,newaxis]
        # more shortcuts
        variance = self._variance
        tmp = self._tmp
        min_var = self._min_var
        # filter buffer for every radius to calc variance
        zz = self._xform(self.z)
        for n, (ra, ri) in enumerate(zip(self.ra, self.ri)):
            self._filter(zz, ra, output=variance[n,...])
            self._filter(zz, ri, output=tmp)
            variance[n,...] -= tmp
        # calc minimum variance
        min_var = argmin(variance ** 2, axis=0)
        # update greyscale buffer z
#        self.z += choose(min_var, self.dt[_XY] * sign(variance)) * speed
        self.z += speed * .01 * sum(sign(variance)/variance**2,axis=0) / sum(1/variance**2,axis=0)
        # normalize
        mi, ma = self.z.min(), self.z.max()
        self.z = (self.z - mi) / (ma - mi)
        # update colour buffer c
        # self.c = (1-self._dc) * self.c + self._dc * choose(min_var[...,newaxis], self.pal[:,newaxis,newaxis,:])
    
    def rgb_image(self):
        """Return RGB image from greyscale buffer z combined with colour buffer c."""
        z3 = self.z[:,:,newaxis]
        return z3 * self.c
    
    def run(self, n=1, speed=1.0, rnd=0, filename=None, start_frame=0, verbose=True, crop=None):
        """Advance the multiscale Turing pattern by n timesteps.

        Keyword arguments:
        n -- number of timesteps
        speed -- dt multiplier, default 1.0
        rnd -- noise to add at each step, range 0..1, default 0
        filename -- filename pattern to save intermediate frames
            e.g. 'MSTP%04d.png'
        start_frame -- number of the first frame with respect to the
            filename. Useful for continuing an interrupted run sequence.
        verbose -- if True (default), show countdown during rendering
        crop -- rect to crop the saved image, default None, meaning no cropping.

        """
        if verbose and filename:
            print 'rendering %s frames as %s ... %s' % (n, (filename % start_frame), (filename % (start_frame + n - 1)))
        for k in xrange(n):
            self.z += rnd * rand(*self.z.shape)
            self.step(speed=speed)
            if filename:
                out = self.rgb_image()
                if crop:
                    out = out[crop[0]:crop[1],crop[2]:crop[3],...]
                imsave(filename % (k + start_frame), out)
            if verbose:
                print n - k,
                sys.stdout.flush()

    def imshow(self):
        """pylab.imshow of self.rgb_image with display area maximized by removing margins and tickmarks."""
        axes([0, 0, 1, 1], xticks=[], yticks=[])
        imshow(self.rgb_image())

    def __repr__(self):
        """Optimized but complete representation self (excludes buffers)."""
        def r(a):
            rr = re.sub('\s','',repr(a))
            rr = re.sub('\.(?=\D)','',rr)
            rr = re.sub('(?<=\D)0\.','.',rr)
            if rr.startswith('array(['):
                rr = rr[6:-1]
            return rr
        return 'MSTP(shape=%s,\n    ra=%s,\n    ri=%s,\n    dt=%s,\n    pal=%s)' % (
                r(self.z.shape), r(self.ra), r(self.ri), r(self.dt), r(self.pal))
                
    def __call__(self):
        display(IM(self.rgb_image()))

def shade(N=10, *colors):    
    """Interpolates colors to construct a smooth palette.

    N -- number of entries in palette
    colors -- the color points, each must have 3 elements
    """
    return interp1d(linspace(0, 1, len(colors)), colors, axis=0)(linspace(0, 1, N)).round(2)
#    cc = array(colors)
#    
#    x = frange(0, 1, npts=N)
#         linspace(0, 1, len(colors))
#    xp = frange(0, 1, npts=len(colors))
#    return array([
#        interp(x,xp,cc[:,0]),
#        interp(x,xp,cc[:,1]),
#        interp(x,xp,cc[:,2])]).transpose().round(2)

def rseq(start=0.0, stop=1.0, N=10, randomness=0.5):
    """Generate monotonically increasing random sequence.

    start, stop -- range, inclusive. like frange
    N -- number of points
    randomness -- 0 is completely regular, equivalent to 
        frange(start, stop, npts=N). 1 is completely random (but
        still monotonically increasing).
    
    Examples:
    >>> rseq(1, 4, N=4, randomness=0.1)     #doctest: +SKIP
    array([ 1.04117853,  1.95102362,  2.98915379,  3.95667588])
    >>> rseq(1, 4, N=4, randomness=1)       #doctest: +SKIP
    array([ 1.21479822,  1.37313434,  2.58511447,  3.63824901])
    """

    return (randomness * sort(start + (stop - start) * rand(N)) 
        + (1 - randomness) * frange(start, stop, npts=N))   

def run_variations(filename, N=50, N_scales=7, shape=(512, 512), 
        steps=[(10,4.0), (20,1.0)], display_inline=False, min_radius=1.5,
        max_radius=90.0):
    """Generates random MSTP parameter configurations and renders them.

    filename -- path and filename for pylab.imsave to write output to.
        must include '%s' or '%03d' or similar for numbering.
    N -- number of variations to generate.
    N_scales -- number of scales, length of ra, ri, dt and pal MSTP 
        parameters.
    shape -- dimensions of MSTP and output images
    steps -- list of (N_steps, relative_speed) tuples, defaults to 10
        steps at 4x speed followed by 20 steps at 1x speed.
    min_radius, max_radius -- lower and upper limits of activator radius 
    display_inline -- call IPython QTConsole function display to 
        display intermediate results inline.

    TODO: more parameters to adjust other particulars of random 
        configuration generation. 

    """

    colors = [[1, 0, 0], [0, 1, 0], [0, 0, 0.9], [1, 1, 0], 
              [1, 0, 1], [0, 1, 1], [1, 1, 1], [1, 0.6, 0]]
    perm = np.random.permutation
    for n in range(N):
        # generate palette
        pal = perm(colors)[:4]
        pal = perm(r_[pal, pal])
        pal = (pal * .6 + roll(pal, 4, axis=0) * 0.4)[:4]
        pal = shade(N_scales, *pal)
        # generate parameters and MSTP object
        ra = exp(rseq(log(min_radius), log(max_radius), N=N_scales, randomness=0.9)).round(2)
        ri = ((1.333 + rand(N_scales)) * ra).round(2)
        dt = (.01 * frange(1, N_scales) ** 0.8).round(3)
        wt = (1.33 + arctan(5 * (rand(N_scales) - .5))).round(2)
        m = MSTP(shape, ra=ra, ri=ri, dt=dt, wt=wt, pal=pal)
        
        print '\n-------- rendering image', filename % n
        print m
        # display(HTML(' '.join('<span style="background:rgb(%d,%d,%d)">(o_O)</span> ' % tuple(255*k) for k in array(pal))))
        for i, (N_steps, speed) in enumerate(steps):
            print 'rendering %s steps at %1.1fx speed' % (N_steps, dt_multiplier)
            sys.stdout.flush()
            m.run(n=N_steps, speed=speed)
            if display_inline:
                display(IM(m.rgb_image()))
            if i < 0:
                first = False
                print 'renoising after iter 1.',
                m.z += 0.25 * rand(*m.z.shape) * m.z.ptp()

        # display(IM(m.rgb_image()))
        imsave(filename % n, m.rgb_image())

class HTML(object):
    """Helper class for IPython QTConsole to display HTML text inline."""
    def __init__(self, html):
        self.html=html
    def _repr_html_(self):
        return self.html
    
class IM(object):
    """Helper class for IPython QTConsole to display arrays inline."""
    def __init__(self, ar):
        self.png_stream = StringIO()
        toimage(ar).save(self.png_stream, format='png')
    def _repr_png_(self):
        return self.png_stream.getvalue()

class PNG(object):
    """Helper class for IPython QTConsole to display PNG files inline."""
    def __init__(self, fn):
        self.fn=fn
    def _repr_png_(self):
        return open(self.fn).read()


