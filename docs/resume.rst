***************
Resume Last Run
***************

**CDNest** automatically saves the sampling status every **num_max_saves/5** steps into a file named **restart_dnest.txt__xxx**,
where xxx represents the number of steps when the status is saved. **CDNest** can resume the last run using these saved files.

In Python runs, one can keep the Python script unchanged, but only replace **run()** with **restart()** as 

.. code-block:: python
    
    # ... (python code)
    # keep the Python script unchanged but 
    # only replace run() with post_run()

    logz = sampler.restart(restart_file="restart_dnest.txt_10000")

    # ... (python code)

where **CDNest** will resume from the status at 10000 steps. 

Note that **CDNest** does not change the already saved output files 
(such as **sample.txt** and **sample_info.txt**) and just append the newly generated parameter sets to these files. 
The user needs to appropriately handle with the cases when the number of samples in **sample.txt** and **sample_info.txt**
is not the same as the number of steps at which **CDNest** resume running (10000 in the above example). 
Usually, one needs to delete some rows in **sample.txt** and **sample_info.txt** to ensure the number of rows exactly same as 
the number of steps at which **CDNest** resume running (10000 in the above sample).

The full arguments of restart function is 

.. code-block:: python
    
    sampler.restart(restart_file="restart_dnest.txt_10000", max_num_saves=15000, max_num_levels = 20)

where **max_num_saves** controls the maximum number of saves and **max_num_levels** controls the maximum number of levels.