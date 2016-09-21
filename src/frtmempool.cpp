// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2015 The Bitcoin Core developers
// Distributed under the MIT software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "frtmempool.h"

#include "clientversion.h"
#include "consensus/consensus.h"
#include "consensus/validation.h"
#include "main.h"
#include "policy/fees.h"
#include "policy/policy.h"
#include "serialize.h"
#include "streams.h"
#include "timedata.h"
#include "util.h"
#include "utilmoneystr.h"
#include "utiltime.h"
#include "version.h"

using namespace std;
//TODO: FRT_WEIGHT, FRT_USAGE_SIZE, FRT_SIZE to be decided
CFrtMemPoolEntry::CFrtMemPoolEntry(const CBlockHeader& _frt, //const CAmount& _nFee,
    int64_t _nTime,
    /*double _entryPriority,*/ unsigned int _entryHeight //,
                                                         // bool poolHasNoInputsOf, CAmount _inChainInputValue,
    /* bool _spendsCoinbase, int64_t _sigOpsCost, LockPoints lp*/) : frt(std::make_shared<CBlockHeader>(_frt)), /*nFee(_nFee),*/ nTime(_nTime), /*entryPriority(_entryPriority),*/ entryHeight(_entryHeight)
//    hadNoDependencies(poolHasNoInputsOf), inChainInputValue(_inChainInputValue),
/*    spendsCoinbase(_spendsCoinbase),*/ //sigOpCost(_sigOpsCost)//, lockPoints(lp)
{
    nFrtWeight = GetFrtSize(); //FRT_WEIGHT;    //GetTransactionWeight(_frt);	//TODO: used in miner.cpp for block weight
                               //    nModSize = _tx.CalculateModifiedSize(GetTxSize());
    nUsageSize = GetFrtSize(); //FRT_USAGESIZE; //RecursiveDynamicUsage(*tx) + memusage::DynamicUsage(tx); //TODO: call func in core_memusage and memusage  -> tx.vin, tx.vout, tx.wit

    /*   nCountWithDescendants = 1;
    nSizeWithDescendants = GetTxSize();
    nModFeesWithDescendants = nFee;
    CAmount nValueIn = _tx.GetValueOut()+nFee;
    assert(inChainInputValue <= nValueIn); */

    //    feeDelta = 0;

    /*   nCountWithAncestors = 1;
    nSizeWithAncestors = GetTxSize();
    nModFeesWithAncestors = nFee;
    nSigOpCostWithAncestors = sigOpCost; */
}

CFrtMemPoolEntry::CFrtMemPoolEntry(const CFrtMemPoolEntry& other)
{
    *this = other;
}
/*
double
CTxMemPoolEntry::GetPriority(unsigned int currentHeight) const
{
    double deltaPriority = ((double)(currentHeight-entryHeight)*inChainInputValue)/nModSize;
    double dResult = entryPriority + deltaPriority;
    if (dResult < 0) // This should only happen if it was called with a height below entry height
        dResult = 0;
    return dResult;
}

void CTxMemPoolEntry::UpdateFeeDelta(int64_t newFeeDelta)
{
    nModFeesWithDescendants += newFeeDelta - feeDelta;
    nModFeesWithAncestors += newFeeDelta - feeDelta;
    feeDelta = newFeeDelta;
}*/
/*
void CFrtMemPoolEntry::UpdateLockPoints(const FruitLockPoints& lp)
{
    lockPoints = lp;
}*/

size_t CFrtMemPoolEntry::GetFrtSize() const
{
    return ::GetSerializeSize(GetFrt(), SER_NETWORK, PROTOCOL_VERSION);
    //GetVirtualTransactionSize(nTxWeight, sigOpCost); //TODO
}

// Update the given tx for any in-mempool descendants.
// Assumes that setMemPoolChildren is correct for the given tx and all
// descendants.
/*
void CTxMemPool::UpdateForDescendants(txiter updateIt, cacheMap &cachedDescendants, const std::set<uint256> &setExclude)
{
    setEntries stageEntries, setAllDescendants;
    stageEntries = GetMemPoolChildren(updateIt);

    while (!stageEntries.empty()) {
        const txiter cit = *stageEntries.begin();
        setAllDescendants.insert(cit);
        stageEntries.erase(cit);
        const setEntries &setChildren = GetMemPoolChildren(cit);
        BOOST_FOREACH(const txiter childEntry, setChildren) {
            cacheMap::iterator cacheIt = cachedDescendants.find(childEntry);
            if (cacheIt != cachedDescendants.end()) {
                // We've already calculated this one, just add the entries for this set
                // but don't traverse again.
                BOOST_FOREACH(const txiter cacheEntry, cacheIt->second) {
                    setAllDescendants.insert(cacheEntry);
                }
            } else if (!setAllDescendants.count(childEntry)) {
                // Schedule for later processing
                stageEntries.insert(childEntry);
            }
        }
    }
    // setAllDescendants now contains all in-mempool descendants of updateIt.
    // Update and add to cached descendant map
    int64_t modifySize = 0;
    CAmount modifyFee = 0;
    int64_t modifyCount = 0;
    BOOST_FOREACH(txiter cit, setAllDescendants) {
        if (!setExclude.count(cit->GetTx().GetHash())) {
            modifySize += cit->GetTxSize();
            modifyFee += cit->GetModifiedFee();
            modifyCount++;
            cachedDescendants[updateIt].insert(cit);
            // Update ancestor state for each descendant
            mapTx.modify(cit, update_ancestor_state(updateIt->GetTxSize(), updateIt->GetModifiedFee(), 1, updateIt->GetSigOpCost()));
        }
    }
    mapTx.modify(updateIt, update_descendant_state(modifySize, modifyFee, modifyCount));
}*/

// vHashesToUpdate is the set of transaction hashes from a disconnected block
// which has been re-added to the mempool.
// for each entry, look for descendants that are outside hashesToUpdate, and
// add fee/size information for such descendants to the parent.
// for each such descendant, also update the ancestor state to include the parent.
/*
void CFrtMemPool::UpdateFruitsFromBlock(const std::vector<uint256> &vHashesToUpdate)
{
    LOCK(cs);
    // For each entry in vHashesToUpdate, store the set of in-mempool, but not
    // in-vHashesToUpdate transactions, so that we don't have to recalculate
    // descendants when we come across a previously seen entry.
//    cacheMap mapMemPoolDescendantsToUpdate;

    // Use a set for lookups into vHashesToUpdate (these entries are already
    // accounted for in the state of their ancestors)
    std::set<uint256> setAlreadyIncluded(vHashesToUpdate.begin(), vHashesToUpdate.end());

    // Iterate in reverse, so that whenever we are looking at at a transaction
    // we are sure that all in-mempool descendants have already been processed.
    // This maximizes the benefit of the descendant cache and guarantees that
    // setMemPoolChildren will be updated, an assumption made in
    // UpdateForDescendants.
    BOOST_REVERSE_FOREACH(const uint256 &hash, vHashesToUpdate) {
        // we cache the in-mempool children to avoid duplicate updates
//        setEntries setChildren;
        // calculate children from mapNextTx
        frtiter it = mapFrt.find(hash);
        if (it == mapFrt.end()) {
            continue;
        }
        auto iter = mapNextTx.lower_bound(COutPoint(hash, 0));
        // First calculate the children, and update setMemPoolChildren to
        // include them, and update their setMemPoolParents to include this tx.
        for (; iter != mapNextTx.end() && iter->first->hash == hash; ++iter) {
            const uint256 &childHash = iter->second->GetHash();
            txiter childIter = mapTx.find(childHash);
            assert(childIter != mapTx.end());
            // We can skip updating entries we've encountered before or that
            // are in the block (which are already accounted for).
            if (setChildren.insert(childIter).second && !setAlreadyIncluded.count(childHash)) {
                UpdateChild(it, childIter, true);
                UpdateParent(childIter, it, true);
            }
        }
        UpdateForDescendants(it, mapMemPoolDescendantsToUpdate, setAlreadyIncluded);
    }
}*/
/*
bool CTxMemPool::CalculateMemPoolAncestors(const CTxMemPoolEntry &entry, setEntries &setAncestors, uint64_t limitAncestorCount, uint64_t limitAncestorSize, uint64_t limitDescendantCount, uint64_t limitDescendantSize, std::string &errString, bool fSearchForParents  = true ) const
{
    setEntries parentHashes;
    const CTransaction &tx = entry.GetTx();

    if (fSearchForParents) {
        // Get parents of this transaction that are in the mempool
        // GetMemPoolParents() is only valid for entries in the mempool, so we
        // iterate mapTx to find parents.
        for (unsigned int i = 0; i < tx.vin.size(); i++) {
            txiter piter = mapTx.find(tx.vin[i].prevout.hash);
            if (piter != mapTx.end()) {
                parentHashes.insert(piter);
                if (parentHashes.size() + 1 > limitAncestorCount) {
                    errString = strprintf("too many unconfirmed parents [limit: %u]", limitAncestorCount);
                    return false;
                }
            }
        }
    } else {
        // If we're not searching for parents, we require this to be an
        // entry in the mempool already.
        txiter it = mapTx.iterator_to(entry);
        parentHashes = GetMemPoolParents(it);
    }

    size_t totalSizeWithAncestors = entry.GetTxSize();

    while (!parentHashes.empty()) {
        txiter stageit = *parentHashes.begin();

        setAncestors.insert(stageit);
        parentHashes.erase(stageit);
        totalSizeWithAncestors += stageit->GetTxSize();

        if (stageit->GetSizeWithDescendants() + entry.GetTxSize() > limitDescendantSize) {
            errString = strprintf("exceeds descendant size limit for tx %s [limit: %u]", stageit->GetTx().GetHash().ToString(), limitDescendantSize);
            return false;
        } else if (stageit->GetCountWithDescendants() + 1 > limitDescendantCount) {
            errString = strprintf("too many descendants for tx %s [limit: %u]", stageit->GetTx().GetHash().ToString(), limitDescendantCount);
            return false;
        } else if (totalSizeWithAncestors > limitAncestorSize) {
            errString = strprintf("exceeds ancestor size limit [limit: %u]", limitAncestorSize);
            return false;
        }

        const setEntries & setMemPoolParents = GetMemPoolParents(stageit);
        BOOST_FOREACH(const txiter &phash, setMemPoolParents) {
            // If this is a new ancestor, add it.
            if (setAncestors.count(phash) == 0) {
                parentHashes.insert(phash);
            }
            if (parentHashes.size() + setAncestors.size() + 1 > limitAncestorCount) {
                errString = strprintf("too many unconfirmed ancestors [limit: %u]", limitAncestorCount);
                return false;
            }
        }
    }

    return true;
}*/
/*
void CTxMemPool::UpdateAncestorsOf(bool add, txiter it, setEntries &setAncestors)
{
    setEntries parentIters = GetMemPoolParents(it);
    // add or remove this tx as a child of each parent
    BOOST_FOREACH(txiter piter, parentIters) {
        UpdateChild(piter, it, add);
    }
    const int64_t updateCount = (add ? 1 : -1);
    const int64_t updateSize = updateCount * it->GetTxSize();
    const CAmount updateFee = updateCount * it->GetModifiedFee();
    BOOST_FOREACH(txiter ancestorIt, setAncestors) {
        mapTx.modify(ancestorIt, update_descendant_state(updateSize, updateFee, updateCount));
    }
}*/
/*
void CTxMemPool::UpdateEntryForAncestors(txiter it, const setEntries &setAncestors)
{
    int64_t updateCount = setAncestors.size();
    int64_t updateSize = 0;
    CAmount updateFee = 0;
    int64_t updateSigOpsCost = 0;
    BOOST_FOREACH(txiter ancestorIt, setAncestors) {
        updateSize += ancestorIt->GetTxSize();
        updateFee += ancestorIt->GetModifiedFee();
        updateSigOpsCost += ancestorIt->GetSigOpCost();
    }
    mapTx.modify(it, update_ancestor_state(updateSize, updateFee, updateCount, updateSigOpsCost));
}*/
/*
void CTxMemPool::UpdateChildrenForRemoval(txiter it)
{
    const setEntries &setMemPoolChildren = GetMemPoolChildren(it);
    BOOST_FOREACH(txiter updateIt, setMemPoolChildren) {
        UpdateParent(updateIt, it, false);
    }
}*/
/*
void CFrtMemPool::UpdateForRemoveFromMempool(const setEntries &entriesToRemove, bool updateDescendants)
{
    // For each entry, walk back all ancestors and decrement size associated with this
    // transaction
    const uint64_t nNoLimit = std::numeric_limits<uint64_t>::max();
    if (updateDescendants) {
        // updateDescendants should be true whenever we're not recursively
        // removing a tx and all its descendants, eg when a transaction is
        // confirmed in a block.
        // Here we only update statistics and not data in mapLinks (which
        // we need to preserve until we're finished with all operations that
        // need to traverse the mempool).
        BOOST_FOREACH(txiter removeIt, entriesToRemove) {
            setEntries setDescendants;
            CalculateDescendants(removeIt, setDescendants);
            setDescendants.erase(removeIt); // don't update state for self
            int64_t modifySize = -((int64_t)removeIt->GetTxSize());
            CAmount modifyFee = -removeIt->GetModifiedFee();
            int modifySigOps = -removeIt->GetSigOpCost();
            BOOST_FOREACH(txiter dit, setDescendants) {
                mapTx.modify(dit, update_ancestor_state(modifySize, modifyFee, -1, modifySigOps));
            }
        }
    }
    BOOST_FOREACH(txiter removeIt, entriesToRemove) {
        setEntries setAncestors;
        const CFrtMemPoolEntry &entry = *removeIt;
        std::string dummy;
        // Since this is a tx that is already in the mempool, we can call CMPA
        // with fSearchForParents = false.  If the mempool is in a consistent
        // state, then using true or false should both be correct, though false
        // should be a bit faster.
        // However, if we happen to be in the middle of processing a reorg, then
        // the mempool can be in an inconsistent state.  In this case, the set
        // of ancestors reachable via mapLinks will be the same as the set of
        // ancestors whose packages include this transaction, because when we
        // add a new transaction to the mempool in addUnchecked(), we assume it
        // has no children, and in the case of a reorg where that assumption is
        // false, the in-mempool children aren't linked to the in-block tx's
        // until UpdateTransactionsFromBlock() is called.
        // So if we're being called during a reorg, ie before
        // UpdateTransactionsFromBlock() has been called, then mapLinks[] will
        // differ from the set of mempool parents we'd calculate by searching,
        // and it's important that we use the mapLinks[] notion of ancestor
        // transactions as the set of things to update for removal.
        //CalculateMemPoolAncestors(entry, setAncestors, nNoLimit, nNoLimit, nNoLimit, nNoLimit, dummy, false);
        // Note that UpdateAncestorsOf severs the child links that point to
        // removeIt in the entries for the parents of removeIt.
//        UpdateAncestorsOf(false, removeIt, setAncestors);
    }
    // After updating all the ancestor sizes, we can now sever the link between each
    // transaction being removed and any mempool children (ie, update setMemPoolParents
    // for each direct child of a transaction being removed).
    BOOST_FOREACH(txiter removeIt, entriesToRemove) {
        UpdateChildrenForRemoval(removeIt);
    }
}*/
/*
void CTxMemPoolEntry::UpdateDescendantState(int64_t modifySize, CAmount modifyFee, int64_t modifyCount)
{
    nSizeWithDescendants += modifySize;
    assert(int64_t(nSizeWithDescendants) > 0);
    nModFeesWithDescendants += modifyFee;
    nCountWithDescendants += modifyCount;
    assert(int64_t(nCountWithDescendants) > 0);
}*/
/*
void CTxMemPoolEntry::UpdateAncestorState(int64_t modifySize, CAmount modifyFee, int64_t modifyCount, int modifySigOps)
{
    nSizeWithAncestors += modifySize;
    assert(int64_t(nSizeWithAncestors) > 0);
    nModFeesWithAncestors += modifyFee;
    nCountWithAncestors += modifyCount;
    assert(int64_t(nCountWithAncestors) > 0);
    nSigOpCostWithAncestors += modifySigOps;
    assert(int(nSigOpCostWithAncestors) >= 0);
}*/

CFrtMemPool::CFrtMemPool() : nFruitsUpdated(0)
{
    _clear(); //lock free clear

    // Sanity checks off by default for performance, because otherwise
    // accepting transactions becomes O(N^2) where N is the number
    // of transactions in the pool
    nCheckFrequency = 0;

    //    minerPolicyEstimator = new CBlockPolicyEstimator(_minReasonableRelayFee);
    //    minReasonableRelayFee = _minReasonableRelayFee;
}

CFrtMemPool::~CFrtMemPool()
{
    //    delete minerPolicyEstimator;
}
/*
void CTxMemPool::pruneSpent(const uint256 &hashTx, CCoins &coins)
{
    LOCK(cs);

    auto it = mapNextTx.lower_bound(COutPoint(hashTx, 0));

    // iterate over all COutPoints in mapNextTx whose hash equals the provided hashTx
    while (it != mapNextTx.end() && it->first->hash == hashTx) {
        coins.Spend(it->first->n); // and remove those outputs from coins
        it++;
    }
}*/

unsigned int CFrtMemPool::GetFruitsUpdated() const
{
    LOCK(cs);
    return nFruitsUpdated;
}

void CFrtMemPool::AddFruitsUpdated(unsigned int n)
{
    LOCK(cs);
    nFruitsUpdated += n;
}

bool CFrtMemPool::addUnchecked(const uint256& hash, const CFrtMemPoolEntry& entry, /*setEntries &setAncestors, bool fCurrentEstimate*/)
{
    // Add to memory pool without checking anything.
    // Used by main.cpp AcceptToMemoryPool(), which DOES do
    // all the appropriate checks.
    LOCK(cs);
    indexed_fruit_set::iterator newit = mapFrt.insert(entry).first;
    //    mapLinks.insert(make_pair(newit, TxLinks()));

    // Update transaction for any feeDelta created by PrioritiseTransaction
    // txTODO: refactor so that the fee delta is calculated before inserting
    // into mapTx.
    /*    std::map<uint256, std::pair<double, CAmount> >::const_iterator pos = mapDeltas.find(hash);
    if (pos != mapDeltas.end()) {
        const std::pair<double, CAmount> &deltas = pos->second;
        if (deltas.second) {
            mapTx.modify(newit, update_fee_delta(deltas.second));
        }
    }*/

    // Update cachedInnerUsage to include contained transaction's usage.
    // (When we update the entry for in-mempool parents, memory usage will be
    // further updated.)
    cachedInnerUsage += entry.DynamicMemoryUsage();

    /*    const CTransaction& tx = newit->GetTx();
    std::set<uint256> setParentTransactions;
    for (unsigned int i = 0; i < tx.vin.size(); i++) {
        mapNextTx.insert(std::make_pair(&tx.vin[i].prevout, &tx));
        setParentTransactions.insert(tx.vin[i].prevout.hash);
    }*/
    // Don't bother worrying about child transactions of this one.
    // Normal case of a new transaction arriving is that there can't be any
    // children, because such children would be orphans.
    // An exception to that is if a transaction enters that used to be in a block.
    // In that case, our disconnect block logic will call UpdateTransactionsFromBlock
    // to clean up the mess we're leaving here.

    // Update ancestors with information about this tx
    /*    BOOST_FOREACH (const uint256 &phash, setParentTransactions) {
        txiter pit = mapTx.find(phash);
        if (pit != mapTx.end()) {
            UpdateParent(newit, pit, true);
        }
    }
    UpdateAncestorsOf(true, newit, setAncestors);
    UpdateEntryForAncestors(newit, setAncestors);*/

    nFruitsUpdated++;
    totalFruitSize += entry.GetFruitSize();
    //    minerPolicyEstimator->processTransaction(entry, fCurrentEstimate);

    vFrtHashes.emplace_back(hash, newit);
    newit->vFrtHashesIdx = vFrtHashes.size() - 1;

    return true;
}

void CFrtMemPool::removeUnchecked(frtiter it)
{
    const uint256 hash = it->GetFrt().GetHash();
    /*    BOOST_FOREACH(const CTxIn& txin, it->GetTx().vin)
        mapNextTx.erase(txin.prevout);*/

    if (vFrtHashes.size() > 1) {
        vFrtHashes[it->vFrtHashesIdx] = std::move(vFrtHashes.back());
        vFrtHashes[it->vFrtHashesIdx].second->vFrtHashesIdx = it->vFrtHashesIdx;
        vFrtHashes.pop_back();
        if (vFrtHashes.size() * 2 < vFrtHashes.capacity())
            vFrtHashes.shrink_to_fit();
    } else
        vFrtHashes.clear();

    totalFrtSize -= it->GetFrtSize();
    cachedInnerUsage -= it->DynamicMemoryUsage();
    //    cachedInnerUsage -= memusage::DynamicUsage(mapLinks[it].parents) + memusage::DynamicUsage(mapLinks[it].children);
    //    mapLinks.erase(it);
    mapFrt.erase(it);
    nFruitsUpdated++;
    //    minerPolicyEstimator->removeTx(hash);
}

// Calculates descendants of entry that are not already in setDescendants, and adds to
// setDescendants. Assumes entryit is already a tx in the mempool and setMemPoolChildren
// is correct for tx and all descendants.
// Also assumes that if an entry is in setDescendants already, then all
// in-mempool descendants of it are already in setDescendants as well, so that we
// can save time by not iterating over those entries.
//the only "Descendant" is entryit itself
void CFrtMemPool::CalculateDescendants(frtiter entryit, setEntries& setDescendants)
{
    setEntries stage;
    if (setDescendants.count(entryit) == 0) {
        stage.insert(entryit);
    }
    // Traverse down the children of entry, only adding children that are not
    // accounted for in setDescendants already (because those children have either
    // already been walked, or will be walked in this iteration).
    while (!stage.empty()) {
        frtiter it = *stage.begin();
        setDescendants.insert(it);
        stage.erase(it);

        /*        const setEntries &setChildren = GetMemPoolChildren(it);
        BOOST_FOREACH(const txiter &childiter, setChildren) {
            if (!setDescendants.count(childiter)) {
                stage.insert(childiter);
            }
        }*/
    }
}
void CFrtMemPool::removeRecursive(const CBlockHeader& origFrt, std::list<CBlockHeader>& removed)
{
    // Remove transaction from memory pool
    {
        LOCK(cs);
        setEntries frtToRemove;
        frtiter origit = mapFrt.find(origFrt.GetHash());
        if (origit != mapFrt.end()) {
            frtToRemove.insert(origit);
        } else {
            // When recursively removing but origTx isn't in the mempool
            // be sure to remove any children that are in the pool. This can
            // happen during chain re-orgs if origTx isn't re-accepted into
            // the mempool for any reason.
            /*            for (unsigned int i = 0; i < origTx.vout.size(); i++) {
                auto it = mapNextTx.find(COutPoint(origTx.GetHash(), i));
                if (it == mapNextTx.end())
                    continue;
                txiter nextit = mapTx.find(it->second->GetHash());
                assert(nextit != mapTx.end());
                txToRemove.insert(nextit);
            }*/
        }
        setEntries setAllRemoves;
        BOOST_FOREACH (frtiter it, frtToRemove) {
            CalculateDescendants(it, setAllRemoves);
        }
        BOOST_FOREACH (frtiter it, setAllRemoves) {
            removed.push_back(it->GetFrt());
        }
        RemoveStaged(setAllRemoves /*, false*/);
    }
}
/*TODO: seem no need for fruit
void CFrtMemPool::removeForReorg(const CCoinsViewCache *pcoins, unsigned int nMemPoolHeight, int flags)
{
    // Remove transactions spending a coinbase which are now immature and no-longer-final transactions
    LOCK(cs);
    list<CBlock> fruitsToRemove;
    for (indexed_fruit_set::const_iterator it = mapFrt.begin(); it != mapFrt.end(); it++) {
        const CBlock& frt = it->GetFrt();
        FruitLockPoints lp = it->GetLockPoints();
        bool validLP =  TestLockPointValidity(&lp);
        if (!CheckFinalTx(tx, flags) || !CheckSequenceLocks(tx, flags, &lp, validLP)) {
            // Note if CheckSequenceLocks fails the LockPoints may still be invalid
            // So it's critical that we remove the tx and not depend on the LockPoints.
            transactionsToRemove.push_back(tx);
        } else if (it->GetSpendsCoinbase()) {
            BOOST_FOREACH(const CTxIn& txin, tx.vin) {
                indexed_transaction_set::const_iterator it2 = mapTx.find(txin.prevout.hash);
                if (it2 != mapTx.end())
                    continue;
                const CCoins *coins = pcoins->AccessCoins(txin.prevout.hash);
		if (nCheckFrequency != 0) assert(coins);
                if (!coins || (coins->IsCoinBase() && ((signed long)nMemPoolHeight) - coins->nHeight < COINBASE_MATURITY)) {
                    transactionsToRemove.push_back(tx);
                    break;
                }
            }
        }
        if (!validLP) {
            mapTx.modify(it, update_lock_points(lp));
        }
    }
    BOOST_FOREACH(const CTransaction& tx, transactionsToRemove) {
        list<CTransaction> removed;
        removeRecursive(tx, removed);
    }
}*/
/*
void CTxMemPool::removeConflicts(const CTransaction &tx, std::list<CTransaction>& removed)
{
    // Remove transactions which depend on inputs of tx, recursively
    LOCK(cs);
    BOOST_FOREACH(const CTxIn &txin, tx.vin) {
        auto it = mapNextTx.find(txin.prevout);
        if (it != mapNextTx.end()) {
            const CTransaction &txConflict = *it->second;
            if (txConflict != tx)
            {
                removeRecursive(txConflict, removed);
                ClearPrioritisation(txConflict.GetHash());
            }
        }
    }
}*/

/**
 * Called when a block is connected. Removes from mempool and updates the miner fee estimator.
 */
void CFrtMemPool::removeForBlock(const std::vector<CBlockHeader>& vfrt /*, unsigned int nBlockHeight,
                                std::list<CBlock>& conflicts, bool fCurrentEstimate*/)
{
    LOCK(cs);
    std::vector<CFrtMemPoolEntry> entries;
    BOOST_FOREACH (const CBlockHeader& frt, vfrt) {
        uint256 hash = frt.GetHash();

        indexed_fruit_set::iterator i = mapFrt.find(hash);
        if (i != mapFrt.end())
            entries.push_back(*i);
    }
    BOOST_FOREACH (const CBlockHeader& frt, vfrt) {
        frtiter it = mapFrt.find(frt.GetHash());
        if (it != mapFrt.end()) {
            setEntries stage;
            stage.insert(it);
            RemoveStaged(stage /*, true*/);
        }
        //        removeConflicts(tx, conflicts);
        //        ClearPrioritisation(tx.GetHash());
    }
    // After the txs in the new block have been removed from the mempool, update policy estimates
    //    minerPolicyEstimator->processBlock(nBlockHeight, entries, fCurrentEstimate);
    //    lastRollingFeeUpdate = GetTime();
    //    blockSinceLastRollingFeeBump = true;
}

void CFrtMemPool::_clear()
{
    //    mapLinks.clear();
    mapFrt.clear();
    //    mapNextTx.clear();
    totalFrtSize = 0;
    cachedInnerUsage = 0;
    //    lastRollingFeeUpdate = GetTime();
    //    blockSinceLastRollingFeeBump = false;
    //    rollingMinimumFeeRate = 0;
    ++nFruitsUpdated;
}

void CFrtMemPool::clear()
{
    LOCK(cs);
    _clear();
}

void CFrtMemPool::check(/*const CCoinsViewCache *pcoins*/) const
{
    if (nCheckFrequency == 0)
        return;

    if (insecure_rand() >= nCheckFrequency)
        return;

    LogPrint("frtmempool", "Checking frtmempool with %u fruits\n" /* and %u inputs\n"*/, (unsigned int)mapFrt.size() /*, (unsigned int)mapNextTx.size()*/);

    uint64_t checkTotal = 0;
    uint64_t innerUsage = 0;

    //    CCoinsViewCache mempoolDuplicate(const_cast<CCoinsViewCache*>(pcoins));
    //    const int64_t nSpendHeight = GetSpendHeight(mempoolDuplicate);

    LOCK(cs);
    list<const CFrtMemPoolEntry*> waitingOnDependants;
    for (indexed_fruit_set::const_iterator it = mapFrt.begin(); it != mapFrt.end(); it++) {
        unsigned int i = 0;
        checkTotal += it->GetFrtSize();
        innerUsage += it->DynamicMemoryUsage(); //TODO: MemoryUsage is not determined yet.
        const CBlockHeader& frt = it->GetFrt();
        //        txlinksMap::const_iterator linksiter = mapLinks.find(it);
        //        assert(linksiter != mapLinks.end());
        //        const TxLinks &links = linksiter->second;
        //        innerUsage += memusage::DynamicUsage(links.parents) + memusage::DynamicUsage(links.children);
        bool fDependsWait = false;
        //        setEntries setParentCheck;
        //        int64_t parentSizes = 0;
        //        int64_t parentSigOpCost = 0;
        /*        BOOST_FOREACH(const CTxIn &txin, tx.vin) {
            // Check that every mempool transaction's inputs refer to available coins, or other mempool tx's.
            indexed_transaction_set::const_iterator it2 = mapTx.find(txin.prevout.hash);
            if (it2 != mapTx.end()) {
                const CTransaction& tx2 = it2->GetTx();
                assert(tx2.vout.size() > txin.prevout.n && !tx2.vout[txin.prevout.n].IsNull());
                fDependsWait = true;
                if (setParentCheck.insert(it2).second) {
                    parentSizes += it2->GetTxSize();
                    parentSigOpCost += it2->GetSigOpCost();
                }
            } else {
                const CCoins* coins = pcoins->AccessCoins(txin.prevout.hash);
                assert(coins && coins->IsAvailable(txin.prevout.n));
            }
            // Check whether its inputs are marked in mapNextTx.
            auto it3 = mapNextTx.find(txin.prevout);
            assert(it3 != mapNextTx.end());
            assert(it3->first == &txin.prevout);
            assert(it3->second == &tx);
            i++;
        }*/
        //        assert(setParentCheck == GetMemPoolParents(it));
        // Verify ancestor state is correct.
        //        setEntries setAncestors;
        //        uint64_t nNoLimit = std::numeric_limits<uint64_t>::max();
        //        std::string dummy;
        //        CalculateMemPoolAncestors(*it, setAncestors, nNoLimit, nNoLimit, nNoLimit, nNoLimit, dummy);
        //        uint64_t nCountCheck = setAncestors.size() + 1;
        //        uint64_t nSizeCheck = it->GetFrtSize();
        //        CAmount nFeesCheck = it->GetModifiedFee();
        //        int64_t nSigOpCheck = it->GetSigOpCost();

        /*        BOOST_FOREACH(txiter ancestorIt, setAncestors) {
            nSizeCheck += ancestorIt->GetTxSize();
            nFeesCheck += ancestorIt->GetModifiedFee();
            nSigOpCheck += ancestorIt->GetSigOpCost();
        }*/

        //        assert(it->GetCountWithAncestors() == nCountCheck);
        //        assert(it->GetSizeWithAncestors() == nSizeCheck);
        //        assert(it->GetSigOpCostWithAncestors() == nSigOpCheck);
        //        assert(it->GetModFeesWithAncestors() == nFeesCheck);

        // Check children against mapNextTx
        //        CTxMemPool::setEntries setChildrenCheck;
        //        auto iter = mapNextTx.lower_bound(COutPoint(it->GetTx().GetHash(), 0));
        //        int64_t childSizes = 0;
        /*        for (; iter != mapNextTx.end() && iter->first->hash == it->GetTx().GetHash(); ++iter) {
            txiter childit = mapTx.find(iter->second->GetHash());
            assert(childit != mapTx.end()); // mapNextTx points to in-mempool transactions
            if (setChildrenCheck.insert(childit).second) {
                childSizes += childit->GetTxSize();
            }
        }*/
        //        assert(setChildrenCheck == GetMemPoolChildren(it));
        // Also check to make sure size is greater than sum with immediate children.
        // just a sanity check, not definitive that this calc is correct...
        //        assert(it->GetSizeWithDescendants() >= childSizes + it->GetTxSize());

        /*        if (fDependsWait)
            waitingOnDependants.push_back(&(*it));
        else {
            CValidationState state;
            bool fCheckResult = tx.IsCoinBase() ||
                Consensus::CheckTxInputs(tx, state, mempoolDuplicate, nSpendHeight);
            assert(fCheckResult);
            UpdateCoins(tx, mempoolDuplicate, 1000000);
        }*/
    }
    unsigned int stepsSinceLastRemove = 0;
    /*    while (!waitingOnDependants.empty()) {
        const CTxMemPoolEntry* entry = waitingOnDependants.front();
        waitingOnDependants.pop_front();
        CValidationState state;
        if (!mempoolDuplicate.HaveInputs(entry->GetTx())) {
            waitingOnDependants.push_back(entry);
            stepsSinceLastRemove++;
            assert(stepsSinceLastRemove < waitingOnDependants.size());
        } else {
            bool fCheckResult = entry->GetTx().IsCoinBase() ||
                Consensus::CheckTxInputs(entry->GetTx(), state, mempoolDuplicate, nSpendHeight);
            assert(fCheckResult);
            UpdateCoins(entry->GetTx(), mempoolDuplicate, 1000000);
            stepsSinceLastRemove = 0;
        }
    }
    for (auto it = mapNextTx.cbegin(); it != mapNextTx.cend(); it++) {
        uint256 hash = it->second->GetHash();
        indexed_transaction_set::const_iterator it2 = mapTx.find(hash);
        const CTransaction& tx = it2->GetTx();
        assert(it2 != mapTx.end());
        assert(&tx == it->second);
    }*/

    assert(totalFrtSize == checkTotal);
    assert(innerUsage == cachedInnerUsage);
}
/*
bool CFrtMemPool::CompareDepthAndScore(const uint256& hasha, const uint256& hashb)
{
    LOCK(cs);
    indexed_fruit_set::const_iterator i = mapFrt.find(hasha);
    if (i == mapFrt.end()) return false;
    indexed_fruit_set::const_iterator j = mapFrt.find(hashb);
    if (j == mapFrt.end()) return true;
    uint64_t counta = i->GetCountWithAncestors();
    uint64_t countb = j->GetCountWithAncestors();
    if (counta == countb) {
        return CompareTxMemPoolEntryByScore()(*i, *j);
    }
    return counta < countb;
}*/
/*
namespace {
class DepthAndScoreComparator
{
public:
    bool operator()(const CTxMemPool::indexed_transaction_set::const_iterator& a, const CTxMemPool::indexed_transaction_set::const_iterator& b)
    {
        uint64_t counta = a->GetCountWithAncestors();
        uint64_t countb = b->GetCountWithAncestors();
        if (counta == countb) {
            return CompareTxMemPoolEntryByScore()(*a, *b);
        }
        return counta < countb;
    }
};
}*/
// no need to sort
std::vector<CFrtMemPool::indexed_fruit_set::const_iterator> CFrtMemPool::GetSortedDepthAndScore() const
{
    std::vector<indexed_fruit_set::const_iterator> iters;
    AssertLockHeld(cs);

    iters.reserve(mapFrt.size());

    for (indexed_fruit_set::iterator mi = mapFrt.begin(); mi != mapFrt.end(); ++mi) {
        iters.push_back(mi);
    }
    //    std::sort(iters.begin(), iters.end(), DepthAndScoreComparator());
    return iters;
}

void CFrtMemPool::queryHashes(vector<uint256>& vfrtid)
{
    LOCK(cs);
    auto iters = GetSortedDepthAndScore();

    vfrtid.clear();
    vfrtid.reserve(mapFrt.size());

    for (auto it : iters) {
        vfrtid.push_back(it->GetFrt().GetHash());
    }
}

std::vector<FrtMempoolInfo> CFrtMemPool::infoAll() const
{
    LOCK(cs);
    auto iters = GetSortedDepthAndScore();

    std::vector<FrtMempoolInfo> ret;
    ret.reserve(mapFrt.size());
    for (auto it : iters) {
        ret.push_back(FrtMempoolInfo{it->GetSharedFrt(), it->GetTime() /*, CFeeRate(it->GetFee(), it->GetTxSize())*/});
    }

    return ret;
}

std::shared_ptr<const CBlockHeader> CFrtMemPool::get(const uint256& hash) const
{
    LOCK(cs);
    indexed_fruit_set::const_iterator i = mapFrt.find(hash);
    if (i == mapFrt.end())
        return nullptr;
    return i->GetSharedFrt();
}

FrtMempoolInfo CFrtMemPool::info(const uint256& hash) const
{
    LOCK(cs);
    indexed_fruit_set::const_iterator i = mapFrt.find(hash);
    if (i == mapFrt.end())
        return FrtMempoolInfo();
    return FrtMempoolInfo{i->GetSharedFrt(), i->GetTime() /*, CFeeRate(i->GetFee(), i->GetTxSize())*/};
}
/*
CFeeRate CTxMemPool::estimateFee(int nBlocks) const
{
    LOCK(cs);
    return minerPolicyEstimator->estimateFee(nBlocks);
}
CFeeRate CTxMemPool::estimateSmartFee(int nBlocks, int *answerFoundAtBlocks) const
{
    LOCK(cs);
    return minerPolicyEstimator->estimateSmartFee(nBlocks, answerFoundAtBlocks, *this);
}
double CTxMemPool::estimatePriority(int nBlocks) const
{
    LOCK(cs);
    return minerPolicyEstimator->estimatePriority(nBlocks);
}
double CTxMemPool::estimateSmartPriority(int nBlocks, int *answerFoundAtBlocks) const
{
    LOCK(cs);
    return minerPolicyEstimator->estimateSmartPriority(nBlocks, answerFoundAtBlocks, *this);
}*/
/*
bool
CTxMemPool::WriteFeeEstimates(CAutoFile& fileout) const
{
    try {
        LOCK(cs);
        fileout << 109900; // version required to read: 0.10.99 or later
        fileout << CLIENT_VERSION; // version that wrote the file
        minerPolicyEstimator->Write(fileout);
    }
    catch (const std::exception&) {
        LogPrintf("CTxMemPool::WriteFeeEstimates(): unable to write policy estimator data (non-fatal)\n");
        return false;
    }
    return true;
}*/
/*
bool
CTxMemPool::ReadFeeEstimates(CAutoFile& filein)
{
    try {
        int nVersionRequired, nVersionThatWrote;
        filein >> nVersionRequired >> nVersionThatWrote;
        if (nVersionRequired > CLIENT_VERSION)
            return error("CTxMemPool::ReadFeeEstimates(): up-version (%d) fee estimate file", nVersionRequired);

        LOCK(cs);
        minerPolicyEstimator->Read(filein);
    }
    catch (const std::exception&) {
        LogPrintf("CTxMemPool::ReadFeeEstimates(): unable to read policy estimator data (non-fatal)\n");
        return false;
    }
    return true;
}*/
/*
void CTxMemPool::PrioritiseTransaction(const uint256 hash, const string strHash, double dPriorityDelta, const CAmount& nFeeDelta)
{
    {
        LOCK(cs);
        std::pair<double, CAmount> &deltas = mapDeltas[hash];
        deltas.first += dPriorityDelta;
        deltas.second += nFeeDelta;
        txiter it = mapTx.find(hash);
        if (it != mapTx.end()) {
            mapTx.modify(it, update_fee_delta(deltas.second));
            // Now update all ancestors' modified fees with descendants
            setEntries setAncestors;
            uint64_t nNoLimit = std::numeric_limits<uint64_t>::max();
            std::string dummy;
            CalculateMemPoolAncestors(*it, setAncestors, nNoLimit, nNoLimit, nNoLimit, nNoLimit, dummy, false);
            BOOST_FOREACH(txiter ancestorIt, setAncestors) {
                mapTx.modify(ancestorIt, update_descendant_state(0, nFeeDelta, 0));
            }
        }
    }
    LogPrintf("PrioritiseTransaction: %s priority += %f, fee += %d\n", strHash, dPriorityDelta, FormatMoney(nFeeDelta));
}*/
/*
void CTxMemPool::ApplyDeltas(const uint256 hash, double &dPriorityDelta, CAmount &nFeeDelta) const
{
    LOCK(cs);
    std::map<uint256, std::pair<double, CAmount> >::const_iterator pos = mapDeltas.find(hash);
    if (pos == mapDeltas.end())
        return;
    const std::pair<double, CAmount> &deltas = pos->second;
    dPriorityDelta += deltas.first;
    nFeeDelta += deltas.second;
}*/
/*
void CTxMemPool::ClearPrioritisation(const uint256 hash)
{
    LOCK(cs);
    mapDeltas.erase(hash);
}*/
/*
bool CTxMemPool::HasNoInputsOf(const CTransaction &tx) const
{
    for (unsigned int i = 0; i < tx.vin.size(); i++)
        if (exists(tx.vin[i].prevout.hash))
            return false;
    return true;
}*/

/*CCoinsViewMemPool::CCoinsViewMemPool(CCoinsView* baseIn, const CTxMemPool& mempoolIn) : CCoinsViewBacked(baseIn), mempool(mempoolIn) { }

bool CCoinsViewMemPool::GetCoins(const uint256 &txid, CCoins &coins) const {
    // If an entry in the mempool exists, always return that one, as it's guaranteed to never
    // conflict with the underlying cache, and it cannot have pruned entries (as it contains full)
    // transactions. First checking the underlying cache risks returning a pruned entry instead.
    shared_ptr<const CTransaction> ptx = mempool.get(txid);
    if (ptx) {
        coins = CCoins(*ptx, MEMPOOL_HEIGHT);
        return true;
    }
    return (base->GetCoins(txid, coins) && !coins.IsPruned());
}*/
/*
bool CCoinsViewMemPool::HaveCoins(const uint256 &txid) const {
    return mempool.exists(txid) || base->HaveCoins(txid);
}*/
//TODO: Need to be decided
size_t CFrtMemPool::DynamicMemoryUsage() const
{
    LOCK(cs);
    // Estimate the overhead of mapTx to be 15 pointers + an allocation, as no exact formula for boost::multi_index_contained is implemented.
    return memusage::MallocUsage(sizeof(CFrtMemPoolEntry) + 15 * sizeof(void*)) * mapFrt.size() + /*memusage::DynamicUsage(mapNextTx) + memusage::DynamicUsage(mapDeltas) + memusage::DynamicUsage(mapLinks) +*/ memusage::DynamicUsage(vFrtHashes) + cachedInnerUsage;
}

void CFrtMemPool::RemoveStaged(setEntries& stage /*, bool updateDescendants*/)
{
    AssertLockHeld(cs);
    //    UpdateForRemoveFromMempool(stage, updateDescendants);
    BOOST_FOREACH (const frtiter& it, stage) {
        removeUnchecked(it);
    }
}

int CFrtMemPool::Expire(int64_t time)
{
    LOCK(cs);
    indexed_fruit_set::index<entry_time>::type::iterator it = mapFrt.get<entry_time>().begin();
    setEntries toremove;
    while (it != mapFrt.get<entry_time>().end() && it->GetTime() < time) {
        toremove.insert(mapFrt.project<0>(it));
        it++;
    }
    setEntries stage;
    BOOST_FOREACH (frtiter removeit, toremove) {
        CalculateDescendants(removeit, stage);
    }
    RemoveStaged(stage /*, false*/);
    return stage.size();
}
/*
bool CFrtMemPool::addUnchecked(const uint256&hash, const CFrtMemPoolEntry &entry, bool fCurrentEstimate)
{
    LOCK(cs);
//    setEntries setAncestors;
//    uint64_t nNoLimit = std::numeric_limits<uint64_t>::max();
//    std::string dummy;
//    CalculateMemPoolAncestors(entry, setAncestors, nNoLimit, nNoLimit, nNoLimit, nNoLimit, dummy);
    return addUnchecked(hash, entry, setAncestors, fCurrentEstimate);
}*/
/*
void CFrtMemPool::UpdateChild(frtiter entry, frtiter child, bool add)
{
    setEntries s;
    if (add && mapLinks[entry].children.insert(child).second) {
        cachedInnerUsage += memusage::IncrementalDynamicUsage(s);
    } else if (!add && mapLinks[entry].children.erase(child)) {
        cachedInnerUsage -= memusage::IncrementalDynamicUsage(s);
    }
}*/
/*
void CTxMemPool::UpdateParent(txiter entry, txiter parent, bool add)
{
    setEntries s;
    if (add && mapLinks[entry].parents.insert(parent).second) {
        cachedInnerUsage += memusage::IncrementalDynamicUsage(s);
    } else if (!add && mapLinks[entry].parents.erase(parent)) {
        cachedInnerUsage -= memusage::IncrementalDynamicUsage(s);
    }
}*/
/*
const CFrtMemPool::setEntries & CFrtMemPool::GetMemPoolParents(frtiter entry) const
{
    assert (entry != mapTx.end());
    txlinksMap::const_iterator it = mapLinks.find(entry);
    assert(it != mapLinks.end());
    return it->second.parents;
}*/
/*
const CTxMemPool::setEntries & CTxMemPool::GetMemPoolChildren(txiter entry) const
{
    assert (entry != mapTx.end());
    txlinksMap::const_iterator it = mapLinks.find(entry);
    assert(it != mapLinks.end());
    return it->second.children;
}*/
/* seems that there is no +oo to prevent exceed size_limit.
CFeeRate CFrtMemPool::GetMinFee(size_t sizelimit) const {
    LOCK(cs);
    if (!blockSinceLastRollingFeeBump || rollingMinimumFeeRate == 0)
        return CFeeRate(rollingMinimumFeeRate);

    int64_t time = GetTime();
    if (time > lastRollingFeeUpdate + 10) {
        double halflife = ROLLING_FEE_HALFLIFE;
        if (DynamicMemoryUsage() < sizelimit / 4)
            halflife /= 4;
        else if (DynamicMemoryUsage() < sizelimit / 2)
            halflife /= 2;

        rollingMinimumFeeRate = rollingMinimumFeeRate / pow(2.0, (time - lastRollingFeeUpdate) / halflife);
        lastRollingFeeUpdate = time;

        if (rollingMinimumFeeRate < minReasonableRelayFee.GetFeePerK() / 2) {
            rollingMinimumFeeRate = 0;
            return CFeeRate(0);
        }
    }
    return std::max(CFeeRate(rollingMinimumFeeRate), minReasonableRelayFee);
}*/
/*
void CTxMemPool::trackPackageRemoved(const CFeeRate& rate) {
    AssertLockHeld(cs);
    if (rate.GetFeePerK() > rollingMinimumFeeRate) {
        rollingMinimumFeeRate = rate.GetFeePerK();
        blockSinceLastRollingFeeBump = false;
    }
}*/
//no priority thus remove the first fruit. TODO: Other order?
void CFrtMemPool::TrimToSize(size_t sizelimit /*, std::vector<uint256>* pvNoSpendsRemaining*/)
{
    LOCK(cs);

    unsigned nFrtnRemoved = 0;
    //    CFeeRate maxFeeRateRemoved(0);
    while (!mapFrt.empty() && DynamicMemoryUsage() > sizelimit) {
        //        indexed_fruit_set::index<descendant_score>::type::iterator it = mapTx.get<descendant_score>().begin();
        indexed_fruit_set::index<entry_time_fruit>::type::iterator it = mapFrt.get<entry_time_fruit>().begin(); // TODO: Currently using entry_time

        // We set the new mempool min fee to the feerate of the removed set, plus the
        // "minimum reasonable fee rate" (ie some value under which we consider txn
        // to have 0 fee). This way, we don't allow txn to enter mempool with feerate
        // equal to txn which were removed with no block in between.
        //        CFeeRate removed(it->GetModFeesWithDescendants(), it->GetSizeWithDescendants());
        //        removed += minReasonableRelayFee;
        //        trackPackageRemoved(removed);
        //        maxFeeRateRemoved = std::max(maxFeeRateRemoved, removed);

        setEntries stage;
        CalculateDescendants(mapTx.project<0>(it), stage);
        nFrtnRemoved += stage.size();

        /*        std::vector<CBlockHeader> frtn;
        if (pvNoSpendsRemaining) {
            frtn.reserve(stage.size());
            BOOST_FOREACH(frtiter it, stage)
                frtn.push_back(it->GetFrt());
        }*/
        RemoveStaged(stage /*, false*/);
        /*        if (pvNoSpendsRemaining) {
            BOOST_FOREACH(const CBlock& frt, frtn) {
                BOOST_FOREACH(const CTxIn& txin, tx.vin) {
                    if (exists(txin.prevout.hash))
                        continue;
                    auto it = mapNextTx.lower_bound(COutPoint(txin.prevout.hash, 0));
                    if (it == mapNextTx.end() || it->first->hash != txin.prevout.hash)
                        pvNoSpendsRemaining->push_back(txin.prevout.hash);
                }
            }
        }*/
    }

    //    if (maxFeeRateRemoved > CFeeRate(0))
    LogPrint("frtmempool", "Removed %u txn\n" /*, rolling minimum fee bumped to %s\n"*/, nTxnRemoved /*, maxFeeRateRemoved.ToString()*/);
}
