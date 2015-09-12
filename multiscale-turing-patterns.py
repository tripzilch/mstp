import scipy.ndimage as ndi
import sys

class MSTP(object):
    def __init__(self,
            ri=[ 96, 32, 12,4.5,1.5],
            ra=[ 64, 24,  9,  3,  1],
            dt=[.04,.03,.03,.02,.02],
            z=[300,300]):
        """Create a new multiscale turing pattern.

        Keyword arguments:
        ri -- list of inhibitor radii for each scale
        ra -- list of activator radii for each scale
        dt -- list of timestep per for each scale
        z -- initial pattern, np.array or 2-tuple to initialize with
            uniform random numbers: np.random.rand(*z)

        """
        self.z = z if len(z) != 2 else np.random.rand(*z)
        N = self.N = len(ri)
        self.dt = np.reshape(dt, [N,1,1])        
        self.radii = unique(concatenate([ri, ra]))
        self.ii = self.radii.searchsorted(ri).reshape([N,1,1])
        self.ia = self.radii.searchsorted(ra).reshape([N,1,1])

    def dz(self):
        zz = array([ndi.gaussian_filter(self.z, r, mode='wrap') 
            for r in self.radii])
        vr = choose(self.ia, zz) - choose(self.ii, zz)
        return choose(argmin(abs(vr), axis=0), self.dt * sign(vr))

    def step(self, n=1, filename=None, cmap='gray', show_progress=False):
        """Advance the multiscale Turing pattern by n timesteps.

        Keyword arguments:
        n -- number of timesteps
        filename -- filename pattern to save intermediate frames
            e.g. 'MSTP%04d.png'
        cmap -- color map to use when saving

        """
        for k in xrange(n):
            self.z = (self.z - self.z.min()) / self.z.ptp()
            self.z += self.dz()
            if filename:
                imsave(filename % k, self.z, cmap=cmap)
            if show_progress:
                print n - k,
                sys.stdout.flush()

    def imshow(self, cmap='gray'):
        axes([0,0,1,1],xticks=[],yticks=[])
        imshow(self.z, cmap=cmap)



# inh=[ndi.gaussian_filter(z,r,mode='wrap') for r in ir]
# act=array(act)
# inh=array(inh)
# act.shape
#?var
# dc=act-inh
# dc.shape
# dt=array(dt)
# dc
# sorted(set(ar+ir))

