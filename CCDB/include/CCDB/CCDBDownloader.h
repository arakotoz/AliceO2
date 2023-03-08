// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <cstdio>
#include <cstdlib>
#include <uv.h>
#include <curl/curl.h>
#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <unordered_map>

#ifndef ALICEO2_CCDBDOWNLOADER_H
#define ALICEO2_CCDBDOWNLOADER_H

using namespace std;

namespace o2
{
namespace ccdb
{

/*
 Some functions below aren't member functions of CCDBDownloader because both curl and libuv require callback functions which have to be either static or non-member.
 Because non-static functions are used in the functions below, they must be non-member.
*/

/**
 * uv_walk callback which is used to close passed handle.
 *
 * @param handle Handle to be closed.
 * @param arg Argument required by callback template. Is not used in this implementation.
 */
void closeHandles(uv_handle_t* handle, void* arg);

/**
 * Called by CURL in order to open a new socket. Newly opened sockets are assigned a timeout timer and added to socketTimerMap.
 *
 * @param clientp Pointer to the CCDBDownloader instance which controls the socket.
 * @param purpose Purpose of opened socket. This parameter is unused but required by the callback template.
 * @param address Structure containing information about family, type and protocol for the socket.
 */
curl_socket_t opensocketCallback(void* clientp, curlsocktype purpose, struct curl_sockaddr* address);

/**
 * Delete the handle.
 *
 * @param handle Handle assigned to this callback.
 */
void onUVClose(uv_handle_t* handle);

class CCDBDownloader
{
 public:
  /**
   * Timer starts for each socket when its respective transfer finishes, and is stopped when another transfer starts for that handle.
   * When the timer runs out it closes the socket. The period for which socket stays open is defined by socketTimoutMS.
   */
  std::unordered_map<curl_socket_t, uv_timer_t*> socketTimerMap;

  /**
   * The UV loop which handles transfers.
   */
  uv_loop_t* loop;

  std::unordered_map<uv_handle_t*, bool> handleMap;
  // ADD COMMENT

  /**
   * Time for which sockets will stay open after last download finishes
   */
  int socketTimoutMS = 4000;

  /**
   * Max number of handles that can be used at the same time
   */
  int maxHandlesInUse = 3;

  // CCDBDownloader(uv_loop_t uv_loop);
  CCDBDownloader(uv_loop_t* uv_loop = nullptr);
  ~CCDBDownloader();

  /**
   * Perform on a single handle in a blocking manner. Has the same effect as curl_easy_perform().
   *
   * @param handle Handle to be performed on. It can be reused or cleaned after perform finishes.
   */
  CURLcode perform(CURL* handle);

  /**
   * Perform on a batch of handles. Callback will be exectuted in it's own thread after all handles finish their transfers.
   *
   * @param handles Handles to be performed on.
   */
  std::vector<CURLcode>* asynchBatchPerformWithCallback(std::vector<CURL*> handles, bool* completionFlag, void (*cbFun)(void*), void* cbData);

  /**
   * Perform on a batch of handles in a blocking manner. Has the same effect as calling curl_easy_perform() on all handles in the vector.
   * @param handleVector Handles to be performed on.
   */
  std::vector<CURLcode> batchBlockingPerform(std::vector<CURL*> handleVector);

  /**
   * Perform on a batch of handles. Completion flag will be set to true when all handles finish their transfers.
   * @param handleVector Handles to be performed on.
   * @param completionFlag Should be set to false before passing it to this function. Will be set to true after all transfers finish.
   */
  std::vector<CURLcode>* batchAsynchPerform(std::vector<CURL*> handleVector, bool* completionFlag);

  /**
   * Limits the number of parallel connections. Should be used only if no transfers are happening.
   */
  void setMaxParallelConnections(int limit);

  /**
   * Limits the time a socket and its connection will be opened after transfer finishes.
   */
  void setSocketTimoutTime(int timoutMS);

 private:
  /**
   * Indicates whether the loop that the downloader is running on has been created by it or provided externally.
   * In case of external loop, the loop will not be closed after downloader is deleted.
   */
  bool externalLoop;

  /**
   * Current amount of handles which are performed on.
   */
  int handlesInUse = 0;

  /**
   * Multi handle which controlls all network flow.
   */
  CURLM* curlMultiHandle = nullptr;

  /**
   * The timeout clock that is be used by CURL.
   */
  uv_timer_t* timeout;

  /**
   * Queue of handles awaiting their transfers to start.
   */
  std::vector<CURL*> handlesToBeAdded;

  /**
   * Lock protecting the handleToBeAdded queue.
   */
  std::mutex handlesQueueLock;

  /**
   * Thread on which the thread with uv_loop runs.
   */
  std::thread* loopThread;

  /**
   * Vector with reference to callback threads with a flag marking whether they finished running.
   */
  std::vector<std::pair<std::thread*, bool*>> threadFlagPairVector;

  /**
   * Flag used to signall the loop to close.
   */
  bool closeLoop = false;

  /**
   * Types of requests.
   */
  enum RequestType {
    BLOCKING,
    ASYNCHRONOUS,
    ASYNCHRONOUS_WITH_CALLBACK
  };

  /**
   * Information about a socket.
   */
  typedef struct curl_context_s {
    uv_poll_t poll_handle;
    curl_socket_t sockfd = -1;
    CCDBDownloader* CD = nullptr;
  } curl_context_t;

  /**
   * Structure used for CURLMOPT_SOCKETDATA, which gives context for handleSocket
   */
  typedef struct DataForSocket {
    CCDBDownloader* CD;
    CURLM* curlm;
  } DataForSocket;

  /**
   * Structure which is stored in a easy_handle. It carries information about the request which the easy_handle is part of.
   * All easy handles coming from one request have an identical PerformData structure.
   */
  typedef struct PerformData {
    std::condition_variable* cv;
    bool* completionFlag;
    CURLcode* codeDestination;
    void (*cbFun)(void*);
    std::thread* cbThread;
    void* cbData;
    size_t* requestsLeft;
    RequestType type;
  } PerformData;

  /**
   * Called by CURL in order to close a socket. It will be called by CURL even if a timout timer closed the socket beforehand.
   *
   * @param clientp Pointer to the CCDBDownloader instance which controls the socket.
   * @param item File descriptor of the socket.
   */
  static void closesocketCallback(void* clientp, curl_socket_t item);

  /**
   *  Is used to react to polling file descriptors in poll_handle.
   *
   * @param handle Handle assigned to this callback.
   * @param status Used to signal errors.
   * @param events Bitmask used to describe events on the socket.
   */
  static void curlPerform(uv_poll_t* handle, int status, int events);

  /**
   * Check if loop was signalled to close. The handle connected with this callbacks is always active as to prevent the uv_loop from stopping.
   *
   * @param handle uv_handle to which this callbacks is assigned
   */
  static void checkStopSignal(uv_timer_t* handle);

  /**
   * Used by CURL to react to action happening on a socket.
   */
  static int handleSocket(CURL* easy, curl_socket_t s, int action, void* userp, void* socketp);

  /**
   * Asynchronously notify the loop to check its CURL handle queue.
   *
   * @param handle Handle which is assigned to this callback.
   */
  static void asyncUVHandleCheckQueue(uv_async_t* handle);

  /**
   * Close socket assigned to the timer handle.
   *
   * @param handle Handle which is assigned to this callback.
   */
  static void closeSocketByTimer(uv_timer_t* handle);

  /**
   * Start new transfers, terminate expired transfers.
   *
   * @param req Handle which is assigned to this callback.
   */
  static void curlTimeout(uv_timer_t* req);

  /**
   * Free curl context assigned to the handle.
   *
   * @param handle Handle assigned to this callback.
   */
  static void curlCloseCB(uv_handle_t* handle);

  /**
   * Close poll handle assigned to the socket contained in the context and free data within the handle.
   *
   * @param context Structure containing information about socket and handle to be closed.
   */
  static void destroyCurlContext(curl_context_t* context);

  /**
   * Connect curl timer with uv timer.
   *
   * @param multi Multi handle for which the timout will be set
   * @param timeout_ms Time until timeout
   * @param userp Pointer to the uv_timer_t handle that is used for timeout.
   */
  static int startTimeout(CURLM* multi, long timeout_ms, void* userp);

  /**
   * Check if any of the callback threads have finished running and approprietly join them.
   */
  void checkForThreadsToJoin();

  /**
   * Create a new multi_handle for the downloader
   */
  void initializeMultiHandle();

  /**
   * Release resources reserver for the transfer, mark transfer as complete, passe the CURLcode to the destination and launche callbacks if it is specified in PerformData.
   *
   * @param handle The easy_handle for which the transfer completed
   * @param curlCode The code produced for the handle by the transfer
   */
  void transferFinished(CURL* handle, CURLcode curlCode);

  /**
   * Check message queue inside curl multi handle.
   */
  void checkMultiInfo();

  /**
   * Set openSocketCallback and closeSocketCallback with appropriate arguments. Stores data inside the CURL handle.
   */
  void setHandleOptions(CURL* handle, PerformData* data);

  /**
   * Create structure holding information about a socket including a poll handle assigned to it
   *
   * @param socketfd File descriptor of socket for which the structure will be created
   */
  curl_context_t* createCurlContext(curl_socket_t sockfd);

  /**
   * Asynchroniously signal the event loop to check for new easy_handles to add to multi handle.
   */
  void makeLoopCheckQueueAsync();

  /**
   * If multi_handles uses less then maximum number of handles then add handles from the queue.
   */
  void checkHandleQueue();

  /**
   * Start the event loop. This function should be ran in the `loopThread`.
   */
  void runLoop();
};

/**
 * Structure assigned  to a uv_timer_t before adding it to socketTimerMap. It stores the information about the socket connected to the timer.
 */
typedef struct DataForClosingSocket {
  CCDBDownloader* CD;
  curl_socket_t socket;
} DataForClosingSocket;

} // namespace ccdb
} // namespace o2

#endif
