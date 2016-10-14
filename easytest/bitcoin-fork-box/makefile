BITCOIND=bitcoind
BITCOINGUI=bitcoin-qt
BITCOINCLI=bitcoin-cli
B1_FLAGS=
B2_FLAGS=
B3_FLAGS=
B1=-datadir=1 $(B1_FLAGS)
B2=-datadir=2 $(B2_FLAGS)
B3=-datadir=3 $(B3_FLAGS)
BLOCKS=1
ADDRESS=
AMOUNT=
ACCOUNT=
HASH=
HEIGHT=
STRING=
BOOL=true

start:
	$(BITCOIND) $(B1) -daemon
#	$(BITCOIND) $(B2) -daemon
#	$(BITCOIND) $(B3) -daemon

start-gui:
	$(BITCOINGUI) $(B1) &
	$(BITCOINGUI) $(B2) &
	$(BITCOINGUI) $(B3) &

addnode1:
	$(BITCOINCLI) $(B1) addnode $(ADDRESS) onetry

addnode2:
	$(BITCOINCLI) $(B2) addnode $(ADDRESS) onetry

addnode3:
	$(BITCOINCLI) $(B3) addnode $(ADDRESS) onetry

getaddednodeinfo:
	$(BITCOINCLI) $(B1) getaddednodeinfo $(BOOL)
	$(BITCOINCLI) $(B2) getaddednodeinfo $(BOOL)
	$(BITCOINCLI) $(B3) getaddednodeinfo $(BOOL)	

generate:
	$(BITCOINCLI) $(B1) generate $(BLOCKS)

generate1:
	$(BITCOINCLI) $(B1) generate $(BLOCKS)

generate2:
	$(BITCOINCLI) $(B2) generate $(BLOCKS)

generate3:
	$(BITCOINCLI) $(B3) generate $(BLOCKS)

getblockcount:
	$(BITCOINCLI) $(B1) getblockcount
	$(BITCOINCLI) $(B2) getblockcount
	$(BITCOINCLI) $(B3) getblockcount

gettransaction:
	$(BITCOINCLI) $(B1) gettransaction $(HASH)
	$(BITCOINCLI) $(B2) gettransaction $(HASH)
	$(BITCOINCLI) $(B3) gettransaction $(HASH)

gettransaction1:
	$(BITCOINCLI) $(B1) gettransaction $(HASH)

gettransaction2:
	$(BITCOINCLI) $(B2) gettransaction $(HASH)

gettransaction3:
	$(BITCOINCLI) $(B3) gettransaction $(HASH)

getblockhash:
	$(BITCOINCLI) $(B1) getblockhash $(HEIGHT)
	$(BITCOINCLI) $(B2) getblockhash $(HEIGHT)
	$(BITCOINCLI) $(B3) getblockhash $(HEIGHT)

getblock:
	$(BITCOINCLI) $(B1) getblock $(HASH)
	$(BITCOINCLI) $(B2) getblock $(HASH)
	$(BITCOINCLI) $(B3) getblock $(HASH)

signrawtransaction:
	$(BITCOINCLI) $(B1) signrawtransaction $(HASH)

signrawtransaction1:
	$(BITCOINCLI) $(B1) signrawtransaction $(HASH)

signrawtransaction2:
	$(BITCOINCLI) $(B2) signrawtransaction $(HASH)

signrawtransaction3:
	$(BITCOINCLI) $(B3) signrawtransaction $(HASH)

decoderawtransaction:
	$(BITCOINCLI) $(B1) decoderawtransaction $(HASH)

decoderawtransaction1:
	$(BITCOINCLI) $(B1) decoderawtransaction $(HASH)

decoderawtransaction2:
	$(BITCOINCLI) $(B2) decoderawtransaction $(HASH)

decoderawtransaction3:
	$(BITCOINCLI) $(B3) decoderawtransaction $(HASH)

sendrawtransaction:
	$(BITCOINCLI) $(B1) sendrawtransaction $(HASH)

sendrawtransaction1:
	$(BITCOINCLI) $(B1) sendrawtransaction $(HASH)

sendrawtransaction2:
	$(BITCOINCLI) $(B2) sendrawtransaction $(HASH)

sendrawtransaction3:
	$(BITCOINCLI) $(B3) sendrawtransaction $(HASH)

getinfo:
	$(BITCOINCLI) $(B1) getinfo
	$(BITCOINCLI) $(B2) getinfo
	$(BITCOINCLI) $(B3) getinfo

getrawmempool:
	$(BITCOINCLI) $(B1) getrawmempool $(BOOL)
	$(BITCOINCLI) $(B2) getrawmempool $(BOOL)
	$(BITCOINCLI) $(B3) getrawmempool $(BOOL)

getrawfrtmempool:
	$(BITCOINCLI) $(B1) getrawfrtmempool $(BOOL)
	$(BITCOINCLI) $(B2) getrawfrtmempool $(BOOL)
	$(BITCOINCLI) $(B3) getrawfrtmempool $(BOOL)

getrawfrtmempool_used:
	$(BITCOINCLI) $(B1) getrawfrtmempool_used $(BOOL)
	$(BITCOINCLI) $(B2) getrawfrtmempool_used $(BOOL)
	$(BITCOINCLI) $(B3) getrawfrtmempool_used $(BOOL)

getpeerinfo:
	$(BITCOINCLI) $(B1) getpeerinfo
	$(BITCOINCLI) $(B2) getpeerinfo
	$(BITCOINCLI) $(B3) getpeerinfo

getnetworkinfo:
	$(BITCOINCLI) $(B1) getnetworkinfo
	$(BITCOINCLI) $(B2) getnetworkinfo
	$(BITCOINCLI) $(B3) getnetworkinfo

getaccountaddress1:
	$(BITCOINCLI) $(B1) getaccountaddress ""

getaccountaddress2:
	$(BITCOINCLI) $(B2) getaccountaddress ""

getaccountaddress3:
	$(BITCOINCLI) $(B3) getaccountaddress ""

sendtoaddress:
	$(BITCOINCLI) $(B1) sendtoaddress $(ADDRESS) $(AMOUNT)

sendtoaddress1:
	$(BITCOINCLI) $(B1) sendtoaddress $(ADDRESS) $(AMOUNT)

sendtoaddress2:
	$(BITCOINCLI) $(B2) sendtoaddress $(ADDRESS) $(AMOUNT)

sendtoaddress3:
	$(BITCOINCLI) $(B3) sendtoaddress $(ADDRESS) $(AMOUNT)

address:
	$(BITCOINCLI) $(B1) getnewaddress $(ACCOUNT)

stop:
	$(BITCOINCLI) $(B1) stop
#	$(BITCOINCLI) $(B2) stop
#	$(BITCOINCLI) $(B3) stop

clean:
	rm -rf 1/regtest
#	rm -rf 2/regtest
#	rm -rf 3/regtest

cleanregtest:
	rm -rf 1/regtest
#	rm -rf 2/regtest
#	rm -rf 3/regtest

cleantestnet:
	rm -rf 1/testnet
#	rm -rf 2/regtest
#	rm -rf 3/regtest
