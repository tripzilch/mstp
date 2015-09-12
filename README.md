# Multiscale Turing Patterns, in colour.

This project is quite old, I did it in 2011, so I might not remember everything. I was originally inspired when I saw some work by [Jonathan McCabe](https://vimeo.com/jonathanmccabe) (except it wasn't on vimeo but high-resolution stills on Flickr, and they were in black and white / greyscale). When I saw the shapes and textures created by the MSTP algorithm, I was very impressed and incredibly curious how they were done. They were quite unlike anything else I had seen, either in demoscene-style plasma effects, texture-feedback effects or reaction-diffusion systems (I had done some research into Gray-Scott systems back in university).

So I looked into it, tried to find out all about it, there was a great blog article that explained pretty much everything (except how to make it fast). TODO: Find the link to that article :) If you go looking for it, their implementation was written in Processing (Java).

# The code

Here follows a summary of the source files in this repository. Please note that I worked on this project quite a few years ago, but thankfully it seems I wrote some quite detailed docstrings. So that's *good*!

On the other hand, I also used some real crazy slice indexing voodoo in order to squeeze out a bit more performance. This was probably misguided/premature in light of the performance gains of FFT convolution. So if you're not too familiar with Numpy and advanced array indexing, it's going to seem a little bit arcane. So that's *bad*! Sorry.

It also seems that I actually commented out some parts of the colouring code. Even in `multiscale-turing-patterns-color.py` ... I *think* I de-re-uncommented the lines back correctly but I'm currently not able to test if it works properly. So please use your common sense and otherwise refer to the explanation of the algorithm below. I added some commentary marked `## 2015 NOTE` in `multiscale-turing-patterns-color.py` that tries to explain what is going on, to the best of my knowledge. So that's *vague*! Ah well.

## multiscale-turing-patterns.py

This was my original implementation of the MSTP algorithm. It produces images in greyscale, and uses a Gaussian kernel convolution filter from `scipy.ndimage`, which is not particularly fast for these purposes (because large radiuses).

## multiscale-turing-patterns-color.py

This implementation has colour, which is explained below. It also approximates the Gaussian kernel filter by stacking 3 boxblur filters, but only for larger radiuses. I suppose this is somewhat faster, but not nearly as fast as the FFT algorithm explained next.

## mstp-fft.py

This is most probably the implementation you actually want to use.

Generates MSTPs in colour, but uses FFT to perform the convolutions.

This is so much faster it's not even funny, and I feel kinda stupid for not trying this method right away. I was kind of skeptical about how fast FFT would be, it felt like overkill in a way. But the important factor here is that the convolution kernels are rather *large*, and that's when the
`O(N log(N))` time complexity of FFT really shines in contrast to the `O(N*M)` or `O(N^2)` complexity of a regular convolution. This was definitely a learning experience for me.

The FFT convolutions also allowed me to use (flat) circular convolution kernels instead of Gaussian ones, which is the way it was implemented in the article. From what I remember, this produces sharper details than Gaussian kernels.

## The rest

Leftover bits from various experiments during development. Probably not very useful, included for completeness and educational value. But really, the code in `mstp-fft.py` is where you want to be.

There's some output images in the `output images` folder. These include various high-resolution renders, as well as some experimental stuff: My face in wobbly Turing-patterns, and some images from when I used the MSTP as a procedural texture generation for a platformer game I was working on at the time.

# Colour

So how did I do the colours? I keep two buffers.

First the regular greyscale MSTP buffer, called `z` in my code, using the original algorithm (or at least, my version of it).

I also keep a separate RGB colour buffer called `c`. In the MSTP algorithm, every timestep, for each pixel, the `z` buffer is updated according to the scale (layer) with the lowest variance. I choose a nice-looking palette of colours, one colour for every scale-level. Depending on which scale is chosen for a particular pixel, the pixel in the RGB colour buffer is mixed slightly towards the corresponding colour. That way, the colour buffer keeps track of the history of what scales were recently chosen for that pixel.

Finally, the `z` buffer and the `c` buffer are simply multiplied together, for the output image.

And that's it!

# License stuff

Short bit about usage rights and license stuff: I selected "No License" for this Github repo, which means the regular copyright laws hold (like they always do per default).

As far as I'm concerned, if you just use the code and the explanations to learn and build your own implementation that's just great and awesome. Please do let me know if you make anything cool! :)

If you use the code to play with and experiment with for personal projects, that's also perfectly cool, even if you're not learning :-P

However if you intend to *include* the code (or my writing, or images) in a commercial project, you should send me message first. I'll probably be perfectly fine with that as well, I just want to know about it *beforehand*, in this scenario.

Don't use the picture of my face for anything. It's my face. Mine.
