"""Wrapper around a subset of the multiprocessing built-in library
to make running multiple processes easy. Load balancing is left to
the caller."""

import multiprocessing

CPUS = multiprocessing.cpu_count()


def parallelize(task, inputs, cpus):
    """Parallelize some task on an input list across the specified
    number of CPU cores.

    :param task: name of the function to run
    :type task: function
    :param inputs: arguments to call `task` with
    :type inputs: list or list[tuple]
    :param cpus: number of processor cores to use
    :type cpus: int
    """
    results = []

    # Don't do any work if there are no inputs
    if len(inputs) == 0:
        return results

    # User requested some number of CPUS - ensure it is sane.
    # Then, reduce to match len(inputs) if fewer than CPUS.
    if not 1 <= cpus <= CPUS:
        cpus = CPUS
    cpus = min([cpus, len(inputs)])

    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task, item))

    # Start working on the jobs
    results = start_processes(tasks, cpus)

    return results


def worker(input_queue, output_queue):
    """Worker function to run jobs from the input queue until a stop
    signal is reached.

    :param input_queue: the tasks to run
    :type input_queue: multiprocessing.Queue[tuple[function, tuple]]
    :param output_queue: where to put results returned by each task
    :type output_queue: multiprocessing.Queue
    """
    for func, args in iter(input_queue.get, 'STOP'):
        output_queue.put(func(*args))


def start_processes(inputs, cpus):
    """Spool processes and use them to run jobs.

    :param inputs: the functions and arguments to run in parallel
    :type inputs: list[tuple[func, tuple]]
    :param cpus: the number of CPU cores to use
    :type cpus: int
    :return: results
    """
    job_q = multiprocessing.Queue()
    result_q = multiprocessing.Queue()

    # Populate the job_queue
    [job_q.put(job) for job in inputs]

    # Set up workers and put a 'STOP' signal at the end of job_q for each
    worker_pool = list()
    for i in range(cpus):
        worker_pool.append(
            multiprocessing.Process(target=worker, args=(job_q, result_q)))
        job_q.put('STOP')

    # Ready... set... go!
    [w.start() for w in worker_pool]

    # Grab results from result_q
    results = []
    for _ in range(len(inputs)):
        result = result_q.get()
        if result is not None:
            results.append(result)

    # Make workers join the main process
    [w.join() for w in worker_pool]

    return results
