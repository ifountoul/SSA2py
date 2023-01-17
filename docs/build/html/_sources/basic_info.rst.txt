.. basic_info:

==========
Basic Info
==========

.. raw:: html

    <div class='Basic Info'>

What is SSA2py?
+++++++++++++++

One vital question after the occurrence of a significant seismic event is the spatiotemporal behavior of the seismic source. SSA2py is an emergent python based software that allows fast analysis of the seismic rupture, making possible the near-realtime identification of the rupture characteristics after a significant seismic event. SSA2py does this by implementing the Source-Scanning Algorithm (SSA), a local backprojection implementation suggested by Honn Kao and Shao-Ju Shan (`2004 <https://academic.oup.com/gji/article/157/2/589/633788?login=true>`_, `2007 <https://academic.oup.com/gji/article/168/3/1011/2042108?login=true>`_).

**History of SSA2py**

SSA2py evolved out of resilience scientific work that the NOA-IG team began years ago. The systematic and successful application of SSA together with the urgent need in the seismological community for rapid, near-realtime information about the seismic source led our team to develop SSA2py. Our goal is to provide a continuously evolved software that could contribute substantially to the understanding of the seismic rupture.

The following works provide more context around the SSA methodology, the challenges that it attempts to address and case studies around the world:

- `The Source-Scanning Algorithm: mapping the distribution of seismic sources in time and space <https://academic.oup.com/gji/article/157/2/589/633788?login=true>`_ by Honn Kao and Shao-Ju Shan (2004)
- `Rapid identification of earthquake rupture plane using Source‐Scanning Algorithm <https://academic.oup.com/gji/article/168/3/1011/2042108?login=true>`_ by Honn Kao and Shao-Ju Shan (2007)
- `High-frequency source imaging of the 2011 October 23 Van (Eastern Turkey) earthquake by backprojection of strong motion waveforms <https://academic.oup.com/gji/article/196/2/1060/577351?login=true>`_ by C. P. Evangelidis and H. Kao (2013)
- `Imaging supershear rupture for the 2014 Mw 6.9 Northern Aegean earthquake by backprojection of strong motion waveforms <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014GL062513>`_ by C. P. Evangelidis (2014) 

What is SSA2py For?
++++++++++++++++++++
SSA2py is designed to do the following:

- Enable near real-time monitoring and alerting using FDSN Web Services.
- In case of seismic event the SSA method is applied.
- Parallelized and adapted code to run in GPU and CPU multiprocessing architectures.
- Highly configurable and flexible for more complex applications.
- Extensive scientific analysis of the SSA results.

.. raw:: html

    </div>
