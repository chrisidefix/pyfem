import profile

def main():
  a = 1

profile.run("main()", "profile.tmp")

import pstats
p = pstats.Stats('profile.tmp')
p.sort_stats('cumulative').print_stats(10)
