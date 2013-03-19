import cProfile
import pstats
import model

cProfile.run("model.main()","modelprof")
#p = pstats.Stats('modelprof')
#p.sort_stats('cumulative').print_stats()
#p.sort_stats('time').print_stats()
