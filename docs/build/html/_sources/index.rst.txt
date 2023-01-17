SSA2py Documentation (|release|)
================================

SSA2py is an open-source python project that follows the Source-Scanning Algorithm (SSA), a
back-projection implementation by Honn Kao and Shao-Ju Shan (`2004 <https://academic.oup.com/gji/article/157/2/589/633788?login=true>`_, `2007 <https://academic.oup.com/gji/article/168/3/1011/2042108?login=true>`_). It provides interconnection 
with FDSN Compliant Web Services and it is adapted to run in GPU and CPU multiprocessing architectures.
The aim of SSA2py is to provide rapid and accurate calculations of SSA method in near-realtime conditions.  

.. image:: logo.jpg
  :scale: 30 %


The source code of the project is hosted on `github <https://github.com/ifountoul/SSA2py>`_.

.. raw:: html
    
    <div class="startpage">

Getting Started
+++++++++++++++
     **Basic Info**
         `Few words about SSA2py <basic_info.html>`_ 
     **Methodology**
         `How the program works? <methodology.html>`_
     **Installation**
         `Installation guidelines <installation.html>`_

Start Using SSA2py
++++++++++++++++++
     **Configuration**
         `Details about the configuration of SSA2py <configuration.html>`_
     **Usage**
         `How to use it? <usage.html>`_
     **Case Studies**
         `Examples <case_studies.html>`_

More Info
+++++++++
    **Citation**
         If you used SSA2py in your study, please cite the following **conference presentation**:

         Fountoulakis, I. and Evangelidis, C.: SSA2py: A seismic source imaging tool in Python based on the Source-Scanning Algorithm,
         EGU General Assembly 2022, Vienna, Austria, 23–27 May 2022, EGU22-3800, `<https://doi.org/10.5194/egusphere-egu22-3800>`_, 2022. 

         Hopefully we will have a journal publication soon.

    **Studies Used SSA2py**
        `Detailed list with studies that used SSA2py. <publ_list.html>`_
    **Discussion**
        Concering support, questions, bug fixes and more please submit your post at the Github Discussion board:

        `SSA2py Discussion <https://github.com/ifountoul/SSA2py/discussions>`_

    **Licence**
        Unless defined otherwise in the files, GPLv3 is set. In any case, it's open-source and free to use.

    **Funding**
        - This project is implemented in the `Institute of Geodynamics of the National Observatory of Athens (NOA) <https://bbnet.gein.noa.gr/HL/>`_.
        - This research is financed by the `Hellenic Foundation for Research and Innovation (H.F.R.I.) <https://www.elidek.gr/en/homepage/>`_ under the 
          "First Call for H.F.R.I. Research Projects to support Faculty members and Researchers and the procurement of high-cost research equipment grant" (SIREN project).


    **Contact**
        You can contact us for any question concerning the program:

        - Ioannis Fountoulakis: ifountoul@noa.gr
        - Christos Evangelidis: cevan@noa.gr
 
    **News**
        `News about SSA2py! <news.html>`_
 
 
.. raw:: html

    </div>

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 3
    :caption: Getting Started

    basic_info
    methodology
    installation

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 3
    :caption: Start Using SSA2py

    configuration
    usage
    applications


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 3
    :caption: More Info

    publ_list
    news

