import heapq


class QueueItem:
    def __init__(self, name, priority):
        self.name = name
        self._priority = priority
        self.queue = None

    @property
    def priority(self):
        return self._priority

    @priority.setter
    def priority(self, priority):
        self._priority = priority
        self.queue.prioritize_queue()

    def set_queue(self, queue):
        self.queue = queue

    def __int__(self):
        return self.priority

    def __lt__(self, other):
        return self._priority < other.priority

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class Queue:
    def __init__(self):
        self.queue = []

    def enqueue(self, item: QueueItem):
        heapq.heappush(self.queue, item)
        item.set_queue(self)
        return self

    def prioritize_queue(self):
        heapq.heapify(self.queue)

    def dequeue(self):
        return heapq.heappop(self.queue)

    def __iter__(self):
        return iter(self.queue)

    def __len__(self):
        return len(self.queue)

    def __str__(self):
        return str(self.queue)
