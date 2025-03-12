CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -z noexecstack -march=native -mtune=native \
  -mavx2 -mpopcnt -maes -mbmi2
NISTFLAGS += -Wno-unused-result -O3 -z noexecstack  -mavx2 -mpopcnt -maes -mbmi2 \
  -march=native -mtune=native
SOURCES = sign.c packing.c poly.c reduce.c rounding.c aes.c
HEADERS = config.h params.h api.h sign.h packing.h poly.h \
  reduce.h rounding.h symmetric.h randombytes.h aes.h \
  crypto_declassify.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h
AES_SOURCES = $(SOURCES) fips202.c symmetric-aes.c
AES_HEADERS = $(HEADERS) fips202.h

.PHONY: all speed shared clean

all: \
  fixtest \
  test \
  speed \
  PQCgenKAT_sign1 \
  PQCgenKAT_sign3 \
  PQCgenKAT_sign5 \
  PQCgenKAT_sign5plus 

fixtest: \
  test/test_NIMS1_fix \
  test/test_NIMS3_fix \
  test/test_NIMS5_fix \
  test/test_NIMS5plus_fix 

packing: \
  test/test_packing1 \
  test/test_packing3 \
  test/test_packing5 \
  test/test_packing5plus

test: \
  test/test_NIMS1 \
  test/test_NIMS3 \
  test/test_NIMS5 \
  test/test_NIMS5plus \
  test/test_NIMS1aes \
  test/test_NIMS3aes \
  test/test_NIMS5aes \
  test/test_NIMS5plusaes

speed: \
  test/test_speed1 \
  test/test_speed3 \
  test/test_speed5 \
  test/test_speed5plus \
  test/test_speed1aes \
  test/test_speed3aes \
  test/test_speed5aes \
  test/test_speed5plusaes

timecop: \
  test/test_timecop1 \
  test/test_timecop3 \
  test/test_timecop5 \
  test/test_timecop5plus

shared: \
  libpqcrystals_NIMS1_ref.so \
  libpqcrystals_NIMS3_ref.so \
  libpqcrystals_NIMS5_ref.so \
  libpqcrystals_NIMS5plus_ref.so \
  libpqcrystals_NIMS1aes_ref.so \
  libpqcrystals_NIMS3aes_ref.so \
  libpqcrystals_NIMS5aes_ref.so \
  libpqcrystals_NIMS5plusaes_ref.so \
  libpqcrystals_fips202_ref.so \
  libpqcrystals_aes256ctr_ref.so


test/test_NIMS1: test/test_NIMS.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_NIMS3: test/test_NIMS.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_NIMS5: test/test_NIMS.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)

test/test_NIMS5plus: test/test_NIMS.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_packing1: test/test_packing.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_packing3: test/test_packing.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)

test/test_packing5: test/test_packing.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)

test/test_packing5plus: test/test_packing.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)		

libpqcrystals_fips202_ref.so: fips202.c fips202.h
	$(CC) -shared -fPIC $(CFLAGS) -o $@ $<

libpqcrystals_aes256ctr_ref.so: aes.c aes.h
	$(CC) -shared -fPIC $(CFLAGS) -o $@ $<

libpqcrystals_NIMS1_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $(SOURCES) symmetric-shake.c

libpqcrystals_NIMS3_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $(SOURCES) symmetric-shake.c

libpqcrystals_NIMS5_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $(SOURCES) symmetric-shake.c

libpqcrystals_NIMS5plus_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $(SOURCES) symmetric-shake.c		

libpqcrystals_NIMS1aes_ref.so: $(SOURCES) $(HEADERS) symmetric-aes.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 -DNIMS_TRI_USE_AES \
	   -o $@ $(SOURCES) symmetric-aes.c

libpqcrystals_NIMS3aes_ref.so: $(SOURCES) $(HEADERS) symmetric-aes.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 -DNIMS_TRI_USE_AES \
	   -o $@ $(SOURCES) symmetric-aes.c

libpqcrystals_NIMS5aes_ref.so: $(SOURCES) $(HEADERS) symmetric-aes.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 -DNIMS_TRI_USE_AES \
	   -o $@ $(SOURCES) symmetric-aes.c

libpqcrystals_NIMS5plusaes_ref.so: $(SOURCES) $(HEADERS) symmetric-aes.c
	$(CC) -shared -fPIC $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 -DNIMS_TRI_USE_AES \
	   -o $@ $(SOURCES) symmetric-aes.c		 

test/test_NIMS1aes: test/test_NIMS.c randombytes.c $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 -DNIMS_TRI_USE_AES \
	  -o $@ $< randombytes.c $(AES_SOURCES) 

test/test_NIMS3aes: test/test_NIMS.c randombytes.c $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 -DNIMS_TRI_USE_AES \
	  -o $@ $< randombytes.c $(AES_SOURCES) 

test/test_NIMS5aes: test/test_NIMS.c randombytes.c $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 -DNIMS_TRI_USE_AES \
	  -o $@ $< randombytes.c $(AES_SOURCES)

test/test_NIMS5plusaes: test/test_NIMS.c randombytes.c $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 -DNIMS_TRI_USE_AES \
	  -o $@ $< randombytes.c $(AES_SOURCES)		 

test/test_speed1: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) 

test/test_speed3: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) 

test/test_speed5: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)

test/test_speed5plus: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)		 

test/test_speed1aes: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(AES_SOURCES)  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 -DNIMS_TRI_USE_AES \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(AES_SOURCES) 

test/test_speed3aes: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(AES_SOURCES)  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 -DNIMS_TRI_USE_AES \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(AES_SOURCES) 

test/test_speed5aes: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(AES_SOURCES)  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 -DNIMS_TRI_USE_AES \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(AES_SOURCES)

test/test_speed5plusaes: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(AES_SOURCES)  $(AES_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 -DNIMS_TRI_USE_AES \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(AES_SOURCES)		 

test/test_NIMS1_fix: test/test_NIMS_fixinput.c rng.c rng.h $(KECCAK_SOURCES)  \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) 

test/test_NIMS3_fix: test/test_NIMS_fixinput.c rng.c rng.h $(KECCAK_SOURCES)  \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) 

test/test_NIMS5_fix: test/test_NIMS_fixinput.c rng.c rng.h $(KECCAK_SOURCES)  \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< rng.c $(KECCAK_SOURCES)

test/test_NIMS5plus_fix: test/test_NIMS_fixinput.c rng.c rng.h $(KECCAK_SOURCES)  \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< rng.c $(KECCAK_SOURCES)		 

PQCgenKAT_sign1: PQCgenKAT_sign.c rng.c rng.h $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) 

PQCgenKAT_sign3: PQCgenKAT_sign.c rng.c rng.h $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) 

PQCgenKAT_sign5: PQCgenKAT_sign.c rng.c rng.h $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) 

PQCgenKAT_sign5plus: PQCgenKAT_sign.c rng.c rng.h $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< rng.c $(KECCAK_SOURCES)

PQCgenKAT_sign1aes: PQCgenKAT_sign.c rng.c rng.h $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=1 -DNIMS_TRI_USE_AES \
	  -o $@ $< rng.c $(AES_SOURCES) 

PQCgenKAT_sign3aes: PQCgenKAT_sign.c rng.c rng.h $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=3 -DNIMS_TRI_USE_AES \
	  -o $@ $< rng.c $(AES_SOURCES) 

PQCgenKAT_sign5aes: PQCgenKAT_sign.c rng.c rng.h $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=5 -DNIMS_TRI_USE_AES \
	  -o $@ $< rng.c $(AES_SOURCES)

PQCgenKAT_sign5plusaes: PQCgenKAT_sign.c rng.c rng.h $(AES_SOURCES) \
  $(AES_HEADERS) 
	$(CC) $(NISTFLAGS) -DNIMS_TRI_NTT_MODE=55 -DNIMS_TRI_USE_AES \
	  -o $@ $< rng.c $(AES_SOURCES)		 

test/test_timecop1: test/test_timecop.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=1 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_timecop3: test/test_timecop.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=3 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES) 

test/test_timecop5: test/test_timecop.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=5 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)

test/test_timecop5plus: test/test_timecop.c randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS) 
	$(CC) $(CFLAGS) -DNIMS_TRI_NTT_MODE=55 \
	  -o $@ $< randombytes.c $(KECCAK_SOURCES)		 

clean:
	rm -f *~ test/*~ *.gcno *.gcda *.lcov
	rm -f libpqcrystals_NIMS1_ref.so
	rm -f libpqcrystals_NIMS3_ref.so
	rm -f libpqcrystals_NIMS5_ref.so
	rm -f libpqcrystals_NIMS5plus_ref.so
	rm -f libpqcrystals_NIMS1aes_ref.so
	rm -f libpqcrystals_NIMS3aes_ref.so
	rm -f libpqcrystals_NIMS5aes_ref.so
	rm -f libpqcrystals_NIMS5plusaes_ref.so
	rm -f libpqcrystals_fips202_ref.so
	rm -f libpqcrystals_aes256ctr_ref.so
	rm -f test/test_NIMS1
	rm -f test/test_NIMS3
	rm -f test/test_NIMS5
	rm -f test/test_NIMS5plus
	rm -f test/test_NIMS1aes
	rm -f test/test_NIMS3aes
	rm -f test/test_NIMS5aes
	rm -f test/test_NIMS5plusaes
	rm -f test/test_speed1
	rm -f test/test_speed3
	rm -f test/test_speed5
	rm -f test/test_speed5plus
	rm -f test/test_timecop1
	rm -f test/test_timecop3
	rm -f test/test_timecop5
	rm -f test/test_timecop5plus
	rm -f test/test_speed1aes
	rm -f test/test_speed3aes
	rm -f test/test_speed5aes
	rm -f test/test_speed5plusaes
	rm -f PQCgenKAT_sign1
	rm -f PQCgenKAT_sign3
	rm -f PQCgenKAT_sign5
	rm -f PQCgenKAT_sign5plus
	rm -f PQCgenKAT_sign1aes
	rm -f PQCgenKAT_sign3aes
	rm -f PQCgenKAT_sign5aes
	rm -f PQCgenKAT_sign5plusaes
	rm -f test/test_NIMS1_fix
	rm -f test/test_NIMS3_fix
	rm -f test/test_NIMS5_fix
	rm -f test/test_NIMS5plus_fix
	rm -f test/test_packing1
	rm -f test/test_packing3
	rm -f test/test_packing5
	rm -f test/test_packing5plus
	rm -f *.rsp
	rm -f *.req
