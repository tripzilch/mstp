# Multiscale Turing Patterns, in colour.

Ok so this project is a bit old, I did it in 2011, so I might not remember everything. It was originally inspired when I saw some work by [Jonathan McCabe](https://vimeo.com/jonathanmccabe), except it wasn't on vimeo but high-resolution stills on Flickr, and they were in black and white / greyscale. But when I saw the shapes and textures created by this algorithm, I was very impressed and incredibly curious how they were done. They were quite unlike anything else I had seen, either in demoscene-style plasma effects, texture-feedback effects or reaction-diffusion systems (I had done some research into Gray-Scott systems back in university).

So I looked into it, tried to find out all about it, there was a great blog article that explained pretty much everything (except how to make it fast). TODO: Find the link to that article :) If you go looking for it, their implementation was written in Processing (Java).

# The code

Here follows a quick summary of the source files in this repository. Please note that I worked on this project quite a few years ago, but thankfully it seems I wrote some quite detailed docstrings. So that's good!

On the other hand, I also used some real crazy slice indexing voodoo in order to squeeze out a bit more performance. This was probably misguided/premature in light of the performance gains of FFT convolution. So if you're not too familiar with Numpy and advanced array indexing, it's going to seem a little bit arcane. So that's bad! Sorry.

It also seems that I actually commented out some parts of the colouring code. Even in `multiscale-turing-patterns-color.py` ... I *think* I de-re-uncommented the lines back correctly but I'm currently not able to test if it works properly. So please use your common sense and otherwise refer to the explanation of the algorithm below. I added some commentary marked `## 2015 NOTE` in `multiscale-turing-patterns-color.py` that tries to explain what is going on, to the best of my knowledge. So that's vague! Ah well.

## multiscale-turing-patterns.py

This was my original implementation of the MSTP algorithm. It produces images in greyscale, and uses a Gaussian kernel convolution filter from `scipy.ndimage`, which is not particularly fast for these purposes (because large radiuses).

## multiscale-turing-patterns-color.py

This implementation has colour, which is explained below. It also approximates the Gaussian kernel filter by stacking 3 boxblur filters, but only for larger radiuses. I suppose this is somewhat faster, but not nearly as fast as the FFT algorithm explained next.

## mstp-fft.py

This is most probably the implementation you actually want to use.

Generates MSTPs in colour, but uses FFT to perform the convolutions.

This is so much faster it's not even funny, and I feel kinda stupid for not trying this method right away. I was kind of skeptical about how fast FFT would be, it felt like overkill in a way. But the important factor here is that the convolution kernels are rather *large*, and that's when the
`O(N log(N))` time complexity of FFT really shines in contrast to the `O(N*M)` or `O(N^2)` complexity of a regular convolution.

The FFT convolutions also allowed me to use (flat) circular convolution kernels instead of Gaussian ones, which is the way it was implemented in the article. AFAIK, this produces sharper details than Gaussian kernels.

# Colours

So how did I do the colours?

I keep two buffers. First the regular greyscale MSTP buffer, called `z` in my code, using the original algorithm (or at least, my version of it).

I also keep a separate RGB colour buffer called `c`. In the MSTP algorithm, every timestep, for each pixel, the `z` buffer is updated according to the scale (layer) with the lowest variance. I choose a nice-looking palette of colours, one colour for every scale-level. Depending on which scale is chosen for a particular pixel, the pixel in the RGB colour buffer is mixed slightly towards the corresponding colour. That way, the colour buffer keeps track of the history of what scales were recently chosen for that pixel.

Finally, for the output image, the `z` buffer and the `c` buffer are simply multiplied together.

That's it!













