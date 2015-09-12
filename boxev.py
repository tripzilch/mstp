"""
    Some crazy genetic algorithm / simulated annealing optimization code for trying to synthesize
    combinations of boxblur filters to approximate particular convolution kernels, specifically
    the circular ones used in MSTP.

    This was a crazy experiment idea, somewhere along the road of trying to optimize performance
    of the MSTP convolution filters.

    Long story short, this approach didn't work very well.

    Just use FFT for fast, arbitrary sized convolution filters. It's faster and simpler than anything
    else and I should've tried it first thing, really. Lesson learned!
"""

def sho(y):
    y=y[...,newaxis]
    display(IM(clip(y*[2,1,0],0,1)+clip(-y*[0,1,2],0,1)))

dot = zeros((200,200))
dot[100,100] = 1
tgt = zeros((200,200))
for xx in range(-100,100):
    for yy in range(-100,100):
        tgt[xx+100,yy+100] = 30*30 - ((xx*xx+yy*yy))
tgt = clip(tgt, 0, 1)

dd = zeros((200,200,22))
oo = zeros((200,200))

class FF(object):
    def __init__(self, rr, idx, a, sh):
        self.rr = array(rr)
        self.idx = array(idx)
        self.a = array(a)
        self.sh = array(sh)

    def __call__(self, x):
        global dd,oo
        d = dd[..., :len(self.rr) + 1]
        d[..., 0] = x[...]
        for i, di in enumerate(self.idx):
            nd.uniform_filter(d[..., di], self.rr[i], oo, mode='wrap')
            oo /= oo.max()
            d[..., i + 1] = roll(roll(oo, self.sh[i][0], axis=0), self.sh[i][1], axis=1)
        oo = inner(d, self.a)
        return oo / abs(oo).max()

    def energy(self):
         return mean((self(dot) - tgt)**4.)**.25

    def copy(self):
        return FF(self.rr, self.idx, self.a, self.sh)

    @classmethod
    def random(cls, N=15):
        return cls(
            rr=random_integers(2,50,N),
            idx=[randint(0,i + 1) for i in range(N)],
            a=randn(N + 1),
            sh=random_integers(-5,5,N * 2).reshape(N, 2))

    def dict(self):
        return dict(rr=self.rr, idx=self.idx, a=self.a, sh=self.sh)

    def shake(self, N=2, p=0.5):
        new = self.copy()
        L = len(new.rr)
        for n in range(N):
            if rand() < p:
                i = randint(L+1)
                new.a[i] = .8 * new.a[i] + .2 * randn()
            if rand() < p:
                i = randint(L)
                new.idx[i] = randint(0,i + 1)
            if rand() < p:
                i = randint(L)
                new.rr[i] = (new.rr[i] + randint(9) - 4).clip(2,50)
            if p in [3,8,9]:
                i = randint(L)
                new.sh[i] = (new.sh[i] + randint(5) - 2).clip(-5,5)
        return new

def APF(e, ep, T):
    """Acceptance Probability Function"""
    return 1 if ep < e else exp((e - ep) / T)

fs = FF.random()
fe = fs.energy()
fsbest = fs
febest = fe
print "ANNEALING!!"
for k in range(10000):
    fsnew = fs.shake(N=1)
    fenew = fsnew.energy()
    temp = ((10000.0 - k)/10000)**2 / 200.0
    P = APF(fe, fenew, temp)
    if P > rand():
        if fenew < febest:
            ss = '!!'
        elif fenew < fe:
            ss = '++'
        elif fenew == fe:
            ss = '=='
        else:
            ss = '  '
        fs = fsnew
        fe = fenew
    if fenew < febest:
        fsbest = fsnew
        febest = fenew
        print 'k %04d ------ BEST: %1.4f' % (k, febest)
        sho(hstack([fsbest(dot)[60:140,60:140], (fsbest(dot) - tgt)[60:140,60:140]]))


class BoxGA(object):
    def __init__(self, shape=(200,200), NT=15, pop=200):
        self._shape = shape
        self._NT=NT
        self.d = zeros(shape + (NT + 1,))
        self.x = zeros(shape)
        self.x[shape[0]/2,shape[1]/2] = 1
        self.oo = zeros(shape)
        self.circ = zeros(shape);
        self._pop = 200
        self.circ = clip(self.circ, 0, 1)
        self.reset()

    def reset(self):
        self.tre = [(1.0,dict(rr=self.grr(), idx=self.gid(), a=self.ga(),sh=self.gsh(),corr=1.5+rand())) for i in xrange(self._pop)]
        self.top = []
        self.avg = []

    def ga(self): return randn(self._NT+1)
    def gsh(self): return random_integers(-5,5,self._NT*2).reshape(self._NT,2)
    def grr(self): return random_integers(2,50,self._NT)
    def gid(self): return [randint(0,i + 1) for i in range(self._NT)]

    def display(self):
        row = hstack([self.fff(self.x, **self.tre[i][1])[60:160,60:160] for i in range(9)])
        sho(vstack([row,row - tile(self.circ[60:160,60:160],9)]))

    def info(self):
        print '== GENERATION %03d == %01.4f %01.4f %01.4f %01.4f %01.4f %01.4f %01.4f %01.4f %01.4f == AVG %01.4f' %((len(BG.top),)+tuple(BG.tre[i][0] for i in range(9)) + (mean([tt[0] for tt in BG.tre]),))
        sys.stdout.flush()

    def evolve(self, N=1500):

      def mse(a,b): return mean((a-b)**4.)**.25

      def mingle(u, v):
        L = self._NT
        pi,pr,ps = random_integers(0,L,3)
        pa = random_integers(0,L+1)
        pc = rand()
        return dict(
            a=r_[u['a'][:pa], v['a'][pa:]],
            corr=pc*u['corr']+(1-pc)*v['corr'],
            idx=r_[u['idx'][:pi], v['idx'][pi:]],
            rr=r_[u['rr'][:pr], v['rr'][pr:]],
            sh=r_[u['sh'][:ps], v['sh'][ps:]])


      for gen in range(N):
        self.info()
        if gen % 5 == 0:
            self.display()
        for i in reversed(range(self._pop)):
            if i % 20 == 0:
                print i,
                sys.stdout.flush()
            if i < 10: # keep top 9
                tt = self.tre[i][1]
            elif i < 150: # mutate
                tt = self.tre[randint(0, min(i+1, 100))][1]
                if i % 2 == 0:
                    tt = mingle(tt, self.tre[randint(0, min(i, 100))][1])
                tt = shake(tt, N=4 + i / 25)
            else: # gen new random
                tt = dict(rr=self.grr(), idx=self.gid(), a=self.ga(),sh=self.gsh(),corr=1.5+rand())
                if i < 175:
                    tt = mingle(self.tre[randint(0, i)][1],tt)
            self.tre[i] = (mse(self.fff(self.x, **tt),self.circ),tt)
        self.tre.sort(key=lambda t:t[0])
        self.top.append(self.tre[0][0])
        self.avg.append(mean([tt[0] for tt in self.tre]))

