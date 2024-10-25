*************
Postprocess
*************

**CDNest** output sampling files in real time so that one can diagnostic the sampling results during the running.
In Python runs, this can be done by copying the whole folder to a new elsewhere place, and then calling 
the **post_run()** function, namely, replacing **run()** with **post_run()** in the Python script, 

.. code-block:: python
    
    # ... (python code)
    # keep the Python script unchanged but 
    # only replace run() with post_run()

    logz = sampler.post_run()

    # ... (python code)


If one wants to adjust the posterior temperature, supply a temperature arguement as 

.. code-block:: python
    
    # ... (python code)
    # keep the Python script unchanged but 
    # only replace run() with post_run()

    logz = sampler.post_run(temperature=5)

    # ... (python code)

Here, the posterior temperature is set to 5.