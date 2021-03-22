from __future__ import print_function, division

# Import Python commands
import multiprocessing
import time
import sys


def multi_processer(jobs, cpus=multiprocessing.cpu_count(), info=False, timeout=10, stop_on_error=True, delay=0.3):
    """

    Jobs is expected to be of list/array type and be structured as:
    [  ( functionReference, arguments, key word arguments as dictionary ) , ... ]
    
    stopOnError: [True]  If an error is encountered the processing of all child processes is to be stopped
                 and the error is raised.
                 [False] If an error is encountered while processing, the error will be ignored and the
                 corresponding place in the result vector will be replaced with False.
    
    timeout:     This is a timeout value to avoid dead lock of processes. If the process takes longer than 
                 the prescribed value a timeOutError will be raised and the child process terminated.
    
    info:        Just add some output regarding progress
    
    delay:       Delay between submission of jobs
    """

    # Verify that the number of processes is not more that available    
    if cpus > multiprocessing.cpu_count():
        cpus = multiprocessing.cpu_count()
    if len(jobs) < cpus:
        cpus = len(jobs)

    # Start timer
    start_time = time.time()
    
    try:
        #      Assemble the command for each job
        job_args = []
        for job in jobs:
            function = job[0]
            arguments = job[1] 
            kwarguments = job[2]
            
            if (arguments is None) and (kwarguments is None):
                job_args.append([function, [], {}])
            elif arguments is None:
                job_args.append([function, [], kwarguments])            
            elif kwarguments is None:
                job_args.append([function, arguments, {}])            
            else:                 
                job_args.append([function, arguments, kwarguments])
                
    except:
        print(" ERROR: multiProcessor - The received arguments could not be interpreted")
        print("        The data must be on the following form:")
        print("        [  ( functionReference, arguments, key word arguments as dictionary )\n\t , "
              "( myFun,[x,y,z], {'a':2 ,'b':3} ) \n\t ,  (myFun2, None, None),   ]\n")
        raise

    try:
        # Spawn worker Pool
        worker_pool = multiprocessing.Pool(processes=cpus)

        # Submit jobs to queue
        queue = []
        for job_arg in job_args:
            function = job_arg[0]
            arguments = job_arg[1]
            kwarguments = job_arg[2]
            queue.append(worker_pool.apply_async(function, arguments, kwarguments))
            time.sleep(delay)

        # Get results
        results = []
        
        for i, item in enumerate(queue):
            try:   # noinspection PyBroadException
                results.append(item.get(timeout=timeout))
                if info:
                    print(" Completed %s of %s" % (i+1, len(queue)))
                    sys.stdout.flush()
            except multiprocessing.TimeoutError:
                print("\n ERROR: Timeout\n")
                print("        To avoid dead lock when workers do not operate as intended")
                print("        or an unrecoverable error arises a time out time is set.")
                print(" ")
                print("        The default timeout is 10s. ")
                print(" ")
                print("        This parameter can be manually adjusted as an argument when")
                print("        calling this module, append timeout='number of seconds'")
                print("        to adjust when timeout should occur.")

                if stop_on_error:
                    print("\n\n Terminating child processes")
                    worker_pool.terminate()
                    worker_pool.join()
                    print("-Done\n")
                    raise
                else:
                    results.append(False)
                    continue
            except:
                print("\n\n The following problem was encountered:")
                print(sys.exc_info()[0])  # - Exit type:', sys.exc_info()[0]
                print(sys.exc_info()[1])  # - Exit type:', sys.exc_info()[0]
        
                print("Stop on error: ", stop_on_error)
                if stop_on_error:
                    print("\n\n Terminating child processes")
                    worker_pool.terminate()
                    worker_pool.join()
                    print("-Done\n")
                    raise
                else:
                    results.append(False)
                    continue

        # Close queue and kill workers
        worker_pool.close()
        worker_pool.join()

        end_time = float(round((time.time()-start_time)*10))/10
        if info:
            print(" Program used %s sub processes with a total duration of: %ss" % (cpus, end_time))
        return results
        
    except:
        print("\n\n The following problem was encountered:")
        print(sys.exc_info()[0])  # - Exit type:', sys.exc_info()[0]
        print(sys.exc_info()[1])  # - Exit type:', sys.exc_info()[0]
        raise
