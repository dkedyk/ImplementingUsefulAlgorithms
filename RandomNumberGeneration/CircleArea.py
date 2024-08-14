import math
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
class Circle:
    def __init__(self, p, r):
        self.p = p
        self.r = r
def distance(p1, p2):
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)
def isInCircle(p, circles):
    for c in circles:
        if distance(p, c.p) <= c.r:
            return True
    return False
def findBoundingBox(circles):
    inf = math.inf
    left = inf
    right = -inf
    down = inf
    up = -inf
    for c in circles:
        if c.p.x - c.r < left:
            left = c.p.x - c.r
        if c.p.x + c.r > right:
            right = c.p.x + c.r
        if c.p.y - c.r < down:
            down = c.p.y - c.r
        if c.p.y + c.r > up:
            up = c.p.y + c.r
    return [Point(left, up), Point(right, down)]
import random
def findCoveredArea(circles, nSimulations):
    box = findBoundingBox(circles)
    left = box[0].x
    right = box[1].x
    down = box[1].y
    up = box[0].y
    inCount = 0
    for _ in range(0, nSimulations):
        if isInCircle(Point(random.uniform(left, right),
            random.uniform(down, up)), circles):
            inCount += 1
    return float(inCount)/nSimulations * (up - down) * (right - left)
print(findCoveredArea([Circle(Point(0, 0), 1)], 1000000))#expect Pi
