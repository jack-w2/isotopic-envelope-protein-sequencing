import heapq


class QueueItem:
    def __init__(self, priority, seq, molecular_formula):
        self._priority = priority
        self.seq = seq
        self.molecular_formula = molecular_formula

    @property
    def priority(self):
        return self._priority

    def __int__(self):
        return self.priority

    def __lt__(self, other):
        return self._priority < other.priority

    def __str__(self):
        return f'QueueItem({self.priority}, {self.seq}, {self.molecular_formula}'

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self.molecular_formula)

    def __eq__(self, other):
        return self.molecular_formula == other.molecular_formula


class Queue:
    def __init__(self):
        self.queue = []
        self.visited = set()

    def enqueue(self, item: QueueItem):
        heapq.heappush(self.queue, (item.priority, item))
        item.set_queue(self)
        return self

    def prioritize_queue(self):
        heapq.heapify(self.queue)

    def dequeue(self):
        ret = heapq.heappop(self.queue)[1]
        while ret in self.visited:
            ret = heapq.heappop(self.queue)[1]
        self.visited.add(ret)
        return ret

    def __iter__(self):
        return iter(self.queue)

    def __len__(self):
        return len(self.queue)

    def __str__(self):
        return str(self.queue)
