import imas
import sys
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

#shot = 386; run = 1
#shot = 120; run = 28
shot = 385; run = 30
itm = imas.ids(shot, run)
itm.open()
ids = itm.equilibrium
#twant = 1000.0
twant = 0.0
interpol = 1
ids.getSlice(twant, interpol)

#print ids.time_slice[-1].profiles_1d.psi
#print ids.time_slice[-1].profiles_1d.q
#sys.exit(0)

p = []; pprime = [];
fortdump = open('../transp2imas/fort.313', 'r')
for j in range(1):
  fortdump.readline()
for j in range(100):
  columns = fortdump.readline().strip().split()
  p.append(columns[0]); pprime.append(columns[1]);
fortdump.close()
plt.plot(p, 'b-')
plt.plot(pprime, 'r-')
plt.plot(ids.time_slice[0].profiles_1d.pressure, 'b--')
plt.plot(ids.time_slice[0].profiles_1d.dpressure_dpsi, 'r--')
plt.ylabel('p (blue), p\' (red)')
plt.grid(True)
plt.tight_layout()
plt.savefig('imas-eq-p-prime-initial.png', format = 'png')
plt.show()
plt.close()

F = []; FFprime = [];
fortdump = open('../transp2imas/fort.314', 'r')
for j in range(1):
  fortdump.readline()
for j in range(100):
  columns = fortdump.readline().strip().split()
  F.append(float(columns[0]) - 32.86); FFprime.append(columns[1]);
fortdump.close()
plt.plot(F, 'b-')
plt.plot(FFprime, 'r-')
plt.plot(ids.time_slice[0].profiles_1d.f - 32.86, 'b--')
plt.plot(ids.time_slice[0].profiles_1d.f_df_dpsi, 'r--')
plt.ylabel('F-32.86 (blue), FF\' (red)')
plt.grid(True)
plt.tight_layout()
plt.savefig('imas-eq-F-FFrime-initial.png', format = 'png')
plt.show()
plt.close()

sys.exit(0)

plt.plot(ids.time_slice[-1].profiles_1d.psi, ids.time_slice[-1].profiles_1d.q, 'b')
plt.xlabel('ids.time_slice[-1].profiles_1d.psi')
plt.ylabel('ids.time_slice[-1].profiles_1d.q')
plt.grid(True)
#plt.xlim([0, 1])
plt.tight_layout()
plt.savefig('imas-eq-q.png', format = 'png')
plt.show()
plt.close()

plt.plot(ids.time_slice[-1].profiles_1d.psi, ids.time_slice[-1].profiles_1d.pressure, 'b')
plt.xlabel('ids.time_slice[-1].profiles_1d.psi')
plt.ylabel('ids.time_slice[-1].profiles_1d.pressure')
plt.grid(True)
#plt.xlim([0, 1])
plt.tight_layout()
plt.savefig('imas-eq-p.png', format = 'png')
plt.show()
plt.close()

plt.plot(ids.time_slice[-1].profiles_1d.psi, ids.time_slice[-1].profiles_1d.dpressure_dpsi, 'b')
plt.xlabel('ids.time_slice[-1].profiles_1d.psi')
plt.ylabel('ids.time_slice[-1].profiles_1d.dpressure_dpsi')
plt.grid(True)
#plt.xlim([0, 1])
plt.tight_layout()
plt.savefig('imas-eq-pprime2.png', format = 'png')
plt.show()
plt.close()

plt.plot(ids.time_slice[-1].profiles_1d.psi, ids.time_slice[-1].profiles_1d.f, 'b')
plt.xlabel('ids.time_slice[-1].profiles_1d.psi')
plt.ylabel('ids.time_slice[-1].profiles_1d.f')
plt.grid(True)
#plt.xlim([0, 1])
plt.tight_layout()
plt.savefig('imas-eq-f.png', format = 'png')
plt.show()
plt.close()

plt.plot(ids.time_slice[-1].profiles_1d.psi, ids.time_slice[-1].profiles_1d.f_df_dpsi, 'b')
plt.xlabel('ids.time_slice[-1].profiles_1d.psi')
plt.ylabel('ids.time_slice[-1].profiles_1d.f_df_dpsi')
plt.grid(True)
#plt.xlim([0, 1])
plt.tight_layout()
plt.savefig('imas-eq-FFprime2.png', format = 'png')
plt.show()
plt.close()

print ids.time_slice[-1].profiles_1d.psi[0], ids.time_slice[-1].profiles_1d.psi[1], ids.time_slice[-1].profiles_1d.psi[2]
print ids.time_slice[-1].profiles_1d.pressure[0], ids.time_slice[-1].profiles_1d.pressure[1], ids.time_slice[-1].profiles_1d.pressure[2]
print ids.time_slice[-1].profiles_1d.dpressure_dpsi[0], ids.time_slice[-1].profiles_1d.dpressure_dpsi[1], ids.time_slice[-1].profiles_1d.dpressure_dpsi[2]

print (ids.time_slice[-1].profiles_1d.pressure[2] - ids.time_slice[-1].profiles_1d.pressure[0]) / (ids.time_slice[-1].profiles_1d.psi[2] - ids.time_slice[-1].profiles_1d.psi[0]), ids.time_slice[-1].profiles_1d.dpressure_dpsi[1]

itm.close()
