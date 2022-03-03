"""Wrapper around a subset of the multiprocessing built-in library
to make running multiple processes easy. Load balancing is left to
the caller."""

import multiprocessing

CPUS = multiprocessing.cpu_count()


def parallelize(inputs, cpus, func):
    """Parallelize some task on an input list across the specified
    number of CPU cores.

    :param inputs: arguments to call `func` with
    :type inputs: list
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param func: function to call in parallel
    :type func: function
    :return: results
    """
    # Return early if called with no inputs
    if len(inputs) == 0:
        return list()

    # User requested some number of cores - make sure it's sane
    if cpus < 1 or cpus > CPUS:
        cpus = min([CPUS, len(inputs)])

    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((func, item))

    return start_processes(tasks, cpus)


def worker(input_queue, output_queue):
    """Worker function to run jobs from the input queue until a stop
    signal is reached.

    :param input_queue: thread-safe queue to pull jobs from
    :type input_queue: multiprocessing.Queue
    :param output_queue: thread-safe queue to add results to
    :type output_queue: multiprocessing.Queue
    """
    for func, args in iter(input_queue.get, "STOP"):
        result = func(*args)
        if result:
            output_queue.put(result)


def start_processes(inputs, cpus):
    """Spool up processes and use them to run the jobs.

    :param inputs:
    :type inputs:
    :param cpus:
    :type cpus:
    :return: results
    """
    job_q = multiprocessing.Queue()
    done_q = multiprocessing.Queue()

    [job_q.put(job) for job in inputs]

    worker_pool = list()
    for i in range(cpus):
        worker_pool.append(
            multiprocessing.Process(target=worker, args=(job_q, done_q)))
        job_q.put("STOP")

    [w.start() for w in worker_pool]

    # Grab results
    results = []
    for _ in range(len(inputs)):
        results.append(done_q.get())

    [w.join() for w in worker_pool]

    return results
