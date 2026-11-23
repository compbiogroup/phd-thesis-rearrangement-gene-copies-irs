#!/bin/pypy3
import getopt
import sys
import math

class Genome:
  def __init__(self, _g, _ir):
    # genes
    # self.g = [0] + _g + [len(_g) + 1]
    self.g = _g
    # intergenic regions
    self.ir = [0] + _ir
    # data checking
    if sorted(self.g) != list(range(0, len(self.g))):
      raise ValueError('The genome is not a permutation.')
    if len(self.g) != (len(self.ir)):
      raise Exception('Invalid intergenic regions.')

class Event:
  @classmethod
  def reversal(cls, pi, i, j, x, y):
    # stop case
    if i > j:
      return
    # intergenic parameters restrictions
    if (x < 0) or (y < 0) or (x > pi.ir[i]) or ( y > pi.ir[j + 1]):
      raise Exception('Invalid parameter x or y.')
    # genes
    t = pi.g[i]
    pi.g[i] = pi.g[j]
    pi.g[j] = t
    # intergenic regions
    x_p = (pi.ir[i] - x)
    y_p = (pi.ir[j + 1] - y)
    pi.ir[j + 1] = x_p + y_p
    pi.ir[i] = x + y
    # recursion
    Event.reversal(pi, i + 1, j - 1, 0, pi.ir[j])

  @classmethod
  def transposition(cls, pi, i, j, k, x, y, z):
    # x' and y'
    x_p = pi.ir[i] - x
    y_p = pi.ir[j] - y
    # Transposition simulated using three reversals
    Event.reversal(pi, i, k - 1, x, z)
    Event.reversal(pi, i, i + k - j - 1, x, y_p)
    Event.reversal(pi, i + k - j, k - 1, y, x_p)

class Algorithm:
  def __init__(self, _pi, _iota):
    self.pi = _pi
    self.iota = _iota
    self.breakpoints = []

  def is_breakpoint(self, i):
    if abs(self.pi.g[i] - self.pi.g[i + 1]) != 1:
      return True
    if self.pi.ir[i + 1] != self.iota.ir[max(self.pi.g[i], self.pi.g[i + 1])]:
      return True
    return False

  def is_block(self, a, b):
    if abs(a - b) != 1:
      return False
    p_a = self.pi.g.index(a)
    p_b = self.pi.g.index(b)
    if abs(p_a - p_b) != 1:
      return False
    return self.pi.ir[max(p_a, p_b)] == self.iota.ir[max(a, b)]

  def is_connected(self, i, j):
    # return values:
    # 0 - not connected
    # 1 - connected (case i)
    # 2 - connected (case ii)
    # 3 - connected (case iii)
    # 4 - connected (case iv)
    nucleotides = self.pi.ir[i + 1] + self.pi.ir[j + 1]
    if abs(self.pi.g[i] - self.pi.g[j]) == 1 and (not self.is_block(self.pi.g[i], self.pi.g[j])) and self.iota.ir[max(self.pi.g[i], self.pi.g[j])] <= nucleotides:
      return 1
    if abs(self.pi.g[i + 1] - self.pi.g[j + 1]) == 1 and (not self.is_block(self.pi.g[i + 1], self.pi.g[j + 1])) and self.iota.ir[max(self.pi.g[i + 1], self.pi.g[j + 1])] <= nucleotides:
      return 1
    if abs(self.pi.g[i] - self.pi.g[j + 1]) == 1 and (not self.is_block(self.pi.g[i], self.pi.g[j + 1])) and self.iota.ir[max(self.pi.g[i], self.pi.g[j + 1])] <= nucleotides:
      return 2
    if abs(self.pi.g[i + 1] - self.pi.g[j]) == 1 and (not self.is_block(self.pi.g[i + 1], self.pi.g[j])) and self.iota.ir[max(self.pi.g[i + 1], self.pi.g[j])] <= nucleotides:
      return 3
    if abs(self.pi.g[i] - self.pi.g[i + 1]) == 1 and (not self.is_block(self.pi.g[i], self.pi.g[i + 1])) and self.iota.ir[max(self.pi.g[i], self.pi.g[i + 1])] <= nucleotides:
      return 4
    if abs(self.pi.g[j] - self.pi.g[j + 1]) == 1 and (not self.is_block(self.pi.g[j], self.pi.g[j + 1])) and self.iota.ir[max(self.pi.g[j], self.pi.g[j + 1])] <= nucleotides:
      return 4
    return 0

  def update_breakpoints(self):
    self.breakpoints = []
    for i in range(0, len(self.pi.g) - 1):
      if self.is_breakpoint(i):
        self.breakpoints.append(i)
    return self.breakpoints

# (i, j) / (i + 1, j + 1)
  def case_i(self, i, j):
    if abs(self.pi.g[i] - self.pi.g[j]) == 1:
      goal = self.iota.ir[max(self.pi.g[i], self.pi.g[j])]
      if (self.pi.ir[i + 1] + self.pi.ir[j + 1]) >= goal:
        if self.pi.ir[i + 1] >= goal:
          x = goal
          y = 0
          Event.reversal(self.pi, i + 1, j, x, y)
          return [(i + 1, j, x, y)]
        else:
          x = self.pi.ir[i + 1]
          y = goal - x
          Event.reversal(self.pi, i + 1, j, x, y)
          return [(i + 1, j, x, y)]
    if abs(self.pi.g[i + 1] - self.pi.g[j + 1]) == 1:
      goal = self.iota.ir[max(self.pi.g[i + 1], self.pi.g[j + 1])]
      if (self.pi.ir[i + 1] + self.pi.ir[j + 1]) >= goal:
        if self.pi.ir[j + 1] >= goal:
          x = self.pi.ir[i + 1]
          y = self.pi.ir[j + 1] - goal
          Event.reversal(self.pi, i + 1, j, x, y)
          return [(i + 1, j, x, y)]
        else:
          x = self.pi.ir[i + 1] - (goal - self.pi.ir[j + 1])
          y = 0
          Event.reversal(self.pi, i + 1, j, x, y)
          return [(i + 1, j, x, y)]
    raise Exception('Error - case (i) no operation performed.')

# (i, j + 1)
  def case_ii(self, i, j):
    for k in self.breakpoints:
      if k < i:
        Event.reversal(self.pi, k + 1, i, 0, self.pi.ir[i + 1])
        return [(k + 1, i, 0, self.pi.ir[i + 1])] + self.case_i(k, j)
      if k > j:
        Event.reversal(self.pi, j + 1, k, 0, self.pi.ir[k + 1])
        return [(j + 1, k, 0, self.pi.ir[k + 1])] + self.case_i(i, k)
    raise Exception('Error - case (ii) no breakpoint detected.')

# (i + 1, j)
  def case_iii(self, i, j):
    for k in self.breakpoints:
      if k > i and k < j:
        Event.reversal(self.pi, i + 1, k, 0, self.pi.ir[k + 1])
        return [(i + 1, k, 0, self.pi.ir[k + 1])] + self.case_i(k, j)
    raise Exception('Error - case (iii) no breakpoint detected.')

# (i, i + 1) / (j, j + 1)
  def case_iv(self, i, j):
    Event.reversal(self.pi, i + 1, j, 0, 0)
    return [(i + 1, j, 0, 0)] + self.case_i(i, j)

  def run(self):
    self.update_breakpoints()
    operations = []
    l_b = math.ceil(len(self.breakpoints) / 2.0)
    while(len(self.breakpoints) > 0):
      case = 0
      i = 0
      j = 0
      for a in self.breakpoints:
        for b in self.breakpoints:
          if a < b:
            step = self.is_connected(a, b)
            if step != 0 and (case == 0 or case > step):
              case = step
              i = a
              j = b
      if case == 1:
        operations = operations + self.case_i(i, j)
      if case == 2:
        operations = operations + self.case_ii(i, j)
      if case == 3:
        operations = operations + self.case_iii(i, j)
      if case == 4:
        operations = operations + self.case_iv(i, j)
      self.update_breakpoints()
    # print("{0} {1} {2}".format(int(l_b), len(operations), str(operations).replace(' ', '')))
    print(len(operations))

# if __name__ == "__main__":
#   try :
#     opts, args = getopt.getopt(sys.argv[1:], "hs:",["file="])
#   except getopt.GetoptError:
#     sys.exit()
#   file = None
#   for opt,arg in opts:
#     if(opt == "--file"):
#       file = arg
#   f = open(file, 'r')
#   line = f.readline()
#   line = f.readline()
#   line = f.readline()

#   while line:
#     data = line.split(' ')
#     _pi = eval("{0}".format(data[0]))
#     _pi = [abs(i) for i in _pi]
#     _irpi = eval("{0}".format(data[1]))
#     _iriota = eval("{0}".format(data[2]))

#     pi = Genome(_pi, _irpi)
#     iota = Genome(list(range(1, (len(_pi) + 1))), _iriota)

#     alg = Algorithm(pi, iota)
#     alg.run()
#     line = f.readline()

if __name__ == "__main__":

    _pi = eval("[%s]" % sys.argv[1])
    _pi = [abs(i) for i in _pi]
    _iriota = eval("[%s]" % sys.argv[2])
    _irpi = eval("[%s]" % sys.argv[3])

    pi = Genome(_pi, _irpi)
    iota = Genome(list(range(0, len(_pi))), _iriota)

    alg = Algorithm(pi, iota)
    alg.run()
