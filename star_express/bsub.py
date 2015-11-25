#!/usr/bin/env python3
'''
bsub.py

This script checks a Snakemake job's properties (threads, params) and chooses
an appropriate LSF queue that meets the requirements. It also automatically
chooses the queue that is least busy.

Usage
-----

Add 'threads' and 'params' to your resource-intensive rules:

    rule my_rule:
        input: ...
        output ...
        threads: 4
        params: mem='8000'

Invoke snakemake with the path to this script:

    snakemake --cluster ./bsub.py [other options...]

Note
----

You can use this script with Snakemake on the ERISone cluster at Partners.

For your cluster at your institution, you'll have to modify this script.
'''

import os
import sys

from subprocess import check_output

from snakemake.utils import read_job_properties

def main():
    jobscript = sys.argv[1]

    job_properties = read_job_properties(jobscript)

    # By default, we use 1 thread.
    threads = job_properties.get('threads', 1)

    # We'll leave unspecified the memory and runtime with 0 MB and 0 minutes.
    mem = int(job_properties['params'].get('mem', '0'))
    runtime = int(job_properties['params'].get('runtime', '0'))

    # Let the user specify the queue.
    queue = job_properties['params'].get('queue', None)

    # Otherwise, choose an appropriate queue based on required resources.
    if not queue:
        queue = get_queue(threads, mem, runtime)

    # Submit the job to the queue.
    run_bsub(queue, threads, mem, runtime, jobscript)

def run_bsub(queue, threads, mem, runtime, script):
    cmd = 'bsub -q {q} -n {t}'.format(q=queue, t=threads)
    if mem:
        cmd += ' -R "rusage[mem={m}]"'.format(m=mem)
    if runtime:
        cmd += ' -W {r}'.format(r=runtime)
    cmd += ' {s}'.format(s=script)
    return os.system(cmd)

def get_queue(threads, mem, runtime):
    # All the ERISone queues.
    queues = ['vshort', 'short', 'medium',
              'normal', 'long', 'vlong',
              'big', 'big-multi']
    # Find valid queues for this job's requirements.
    retval = []
    # Only consider 'vshort' if we specify a nonzero runtime.
    if threads == 1 and mem < 1000 and 0 < runtime < 15:
        retval.append('vshort')
    # The other queues are all ok if we leave runtime=0.
    if threads == 1 and mem < 4000 and runtime < 60:
        retval.append('short')
    if threads <= 4 and mem < 8000 and runtime < 60 * 24:
        retval.append('medium')
    if threads <= 6 and mem < 8000 and runtime < 60 * 24 * 3:
        retval.append('normal')
    if threads <= 4 and mem < 8000 and runtime < 60 * 24 * 7:
        retval.append('long')
    if threads <= 4 and mem < 4000 and runtime < 60 * 24 * 7 * 4:
        retval.append('vlong')
    if threads <= 4 and mem > 8000:
        retval.append('big')
    if 4 <= threads <= 16 and mem > 8000:
        retval.append('big-multi')
    # Make sure we have at least one valid queue.
    assert len(retval)
    # Get the number of currently running jobs on each queue.
    lines = check_output('bqueues').split(b'\n')[1:-1]
    lines = [line.decode('utf-8').split() for line in lines]
    njobs = {x[0]: int(x[7]) for x in lines}
    # Among valid queues, choose the one with fewest running jobs.
    return min(retval, key=lambda j: njobs[j])

if __name__ == '__main__':
    main()
