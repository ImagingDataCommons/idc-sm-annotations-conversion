"""Utilities for multiprocessing."""
from multiprocessing import Process, Queue
from typing import Any, Sequence


class _PipelineEnd:
    pass


class Pipeline:

    """A utility class for running sequential operations concurrently.

    A pipeline consists of a sequence of operations, with the output of each
    operation fed to the input of the next.

    The Pipeline runs these processes in parallel on a sequence of data. So
    while operation 2 is processing the first data sample, operation 1 is
    processing the second sample.

    """

    def __init__(
        self,
        operations: Sequence[tuple[type, tuple, dict]],
        maxsize: int = 1,
        same_process: bool = False,
    ):
        """

        Parameters
        ----------
        operations: Sequence[tuple[type, tuple, dict]]
            A list of operation specifications of the form (class, args,
            kwargs).  Each operation is implemented by a callable class (i.e.
            one that implements __class__), and the __call__ method must take a
            single argument. The output of the first class __call__ are fed
            directly to the input of the second, and so on. The args and kwargs
            are passed to the class's __init__ in the sub-process when it is
            first initialized. If any operation returns None, it is not passed
            on to the next operation.
        maxsize: int
            Maximum number of elements that should be stored in each queue at a
            given time.
        same_process: bool
            If True, run the pipeline's operations in the calling thread.

        """
        self._operations = operations
        self._maxsize = maxsize
        self._same_process = same_process

    def __call__(self, data: Sequence[Any]):
        """Run the pipeline on a sequence of data.

        Parameters
        ----------
        data: Sequence[Any]
            Sequence of data to operate on. Each item of this list is passed
            into the first operation in turn.

        """
        if self._same_process:
            processors = [
                cls(*args, **kwargs) for cls, args, kwargs in self._operations
            ]
            for d in data:
                for p in processors:
                    d = p(d)
                    if d is None:
                        break
        else:
            # Processes holds a list of running processes, one for each
            # operation
            processes = []

            # Queues holds a list of queues that connect the processes. There
            # are one fewer queues than processes. The first queue is the
            # output queue of the first output and the input queue of the
            # second operation, and so on.
            queues = []

            # The first operation directly processes data, and therefore has an
            # output queue but no input queue
            src_cls, src_args, src_kwargs = self._operations[0]

            def source_wrapped(data, out_queue):
                obj = src_cls(*src_args, **src_kwargs)
                for d in data:
                    res = obj(d)
                    if res is not None:
                        out_queue.put(res)

                # Use _PipelineEnd to signal that the next process should end
                out_queue.put(_PipelineEnd())

            src_queue = Queue(maxsize=self._maxsize)
            queues.append(src_queue)
            p = Process(target=source_wrapped, args=(data, src_queue))
            p.start()
            processes.append(p)

            # Middle operations (if any) are neither the first or the last, and
            # therefore have both an input queue and an output queue
            for cls, args, kwargs in self._operations[1:-1]:

                def middle_wrapped(in_queue, out_queue):
                    obj = cls(*args, **kwargs)
                    while True:
                        d = in_queue.get()
                        if isinstance(d, _PipelineEnd):
                            break
                        res = obj(d)
                        if res is not None:
                            out_queue.put(res)

                    # Use _PipelineEnd to signal that the next process should
                    # end
                    out_queue.put(_PipelineEnd())

                queues.append(Queue(maxsize=self._maxsize))
                p = Process(
                    target=middle_wrapped,
                    args=(queues[-2], queues[-1])
                )
                p.start()
                processes.append(p)

            # The end operation has no output queue, only an input queue.
            end_cls, end_args, end_kwargs = self._operations[-1]

            def end_wrapped(in_queue):
                obj = end_cls(*end_args, **end_kwargs)
                while True:
                    d = in_queue.get()
                    if isinstance(d, _PipelineEnd):
                        break
                    obj(d)

            p = Process(target=end_wrapped, args=(queues[-1], ))
            p.start()
            processes.append(p)

            for p in processes:
                p.join()
