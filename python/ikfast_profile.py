import pstats

f = open('ikfast_profile.txt', 'w')

p = pstats.Stats('ikfast_profile.out', stream=f)

p.strip_dirs()
p.sort_stats('tottime')
p.print_stats()

p.print_callers()
#p.print_callees()

f.close()


