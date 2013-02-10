#easy_install -U ProDy

import prody as pr

acc = pr.fetchPDB('1ACC')
pacc = pr.parsePDB(acc) #http://www.csb.pitt.edu/ProDy/examples/proteins/parsepdb.html#parsepdb for options

def neighborhood(iterable):
    iterator = iter(iterable)
    prev = None
    item = iterator.next()  # throws StopIteration if empty.
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,None)

for prev,item,next in neighborhood(pacc):
	print prev, item, next
	if ( (prev != None) & (item != None) & (next != None) ):
		print pr.calcAngle(prev, item, next)
		print pr.calcDistance(item, next)