/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Threading;

namespace Numerics.Utilities
{
    /// <summary>
    /// A thread-safe class for reporting the progress of a parallel list of tasks.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     This class was co-developed by Alex Kennedy and Woodrow Fields. Credit to initial implementation and idea goes to Alex.
    ///     Converted to C# and partially documented by Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class SafeProgressReporter
    {

        #region Construction

        /// <summary>
        /// Create a new thread-safe progress reporter.
        /// </summary>
        public SafeProgressReporter()
        {
            _synchronizationContext = SynchronizationContext.Current;
            _invokeProgressHandlers = new SendOrPostCallback(InvokeProgressHandlers);
            _invokeMessageHandlers = new SendOrPostCallback(InvokeMessageHandlers);
        }

        /// <summary>
        /// Create a new thread-safe progress reporter with a specified task name.
        /// </summary>
        /// <param name="taskName">The name of the task.</param>
        public SafeProgressReporter(string taskName)
        {
            TaskName = taskName;
            _synchronizationContext = SynchronizationContext.Current;
            _invokeProgressHandlers = new SendOrPostCallback(InvokeProgressHandlers);
            _invokeMessageHandlers = new SendOrPostCallback(InvokeMessageHandlers);
        }

        #endregion

        #region Members

        private double _previousProgress = -0.0000000000001d;
        private string _previousMessage = "";
        private MessageType _previousMessageType = MessageType.Status;
        private Process _externalProcess;
        private List<SafeProgressReporter> _subProgReporterCollection = new List<SafeProgressReporter>();
        private CancellationTokenSource _cancellationTokenSource = new CancellationTokenSource();
        protected readonly SendOrPostCallback _invokeProgressHandlers;
        protected readonly SendOrPostCallback _invokeMessageHandlers;
        protected SynchronizationContext _synchronizationContext;

        /// <summary>
        /// Returns the most recent progress.
        /// </summary>
        protected double MostRecentProgress => _previousProgress;

        /// <summary>
        /// Returns the most recent message.
        /// </summary>
        public string MostRecentMessage => _previousMessage;

        /// <summary>
        /// Returns the message count.
        /// </summary>
        public int MessageCount { get; private set; }

        /// <summary>
        /// Returns the task name.
        /// </summary>
        public string TaskName { get; private set; }

        /// <summary>
        /// Returns the most recent message type.
        /// </summary>
        public MessageType MostRecentMessageType => _previousMessageType;

        /// <summary>
        /// The external process being executed.
        /// </summary>
        protected Process ExternalProcess => _externalProcess;

        /// <summary>
        /// Determines if cancellation was requested.
        /// </summary>
        public bool IsCancelRequested => _cancellationTokenSource.IsCancellationRequested;

        /// <summary>
        /// Returns a read only list of child progress reporters. 
        /// </summary>
        public ReadOnlyCollection<SafeProgressReporter> ChildReporters
        {
            get { return new ReadOnlyCollection<SafeProgressReporter>(_subProgReporterCollection); }
        }

        /// <summary>
        /// Event is raised when the progress is reported.
        /// </summary>
        public event ProgressReportedEventHandler ProgressReported;
        public delegate void ProgressReportedEventHandler(SafeProgressReporter reporter, double prog, double progDelta);

        /// <summary>
        /// Event is raised when a message is reported.
        /// </summary>
        public event MessageReportedEventHandler MessageReported;
        public delegate void MessageReportedEventHandler(MessageContentStruct msg);

        /// <summary>
        /// Event is raised when the task starts.
        /// </summary>
        public event TaskStartedEventHandler TaskStarted;
        public delegate void TaskStartedEventHandler();

        /// <summary>
        /// Event is raised when the task ended.
        /// </summary>
        public event TaskEndedEventHandler TaskEnded;
        public delegate void TaskEndedEventHandler();

        /// <summary>
        /// Event is raised when a child reporter is created.
        /// </summary>
        public event ChildReporterCreatedEventHandler ChildReporterCreated;
        public delegate void ChildReporterCreatedEventHandler(SafeProgressReporter childReporter);

        /// <summary>
        /// Enumeration of progress reporter message types.
        /// </summary>
        public enum MessageType
        {
            Status,
            Success,
            Warning,
            FatalError
        }

        /// <summary>
        /// Message Content structure.
        /// </summary>
        public struct MessageContentStruct
        {
            public string Message;
            public MessageType MessageType;
            public SafeProgressReporter Reporter;
            public MessageContentStruct(string message, MessageType messageType, SafeProgressReporter reporter)
            {
                Message = message;
                MessageType = messageType;
                Reporter = reporter;
            }
        }

        #endregion

        #region Methods

        /// <summary>
        /// Set synchronization context. 
        /// </summary>
        /// <param name="context">The context.</param>
        protected void SetContext(SynchronizationContext context)
        {
            _synchronizationContext = context;
        }

        /// <summary>
        /// Indicate that the task has started.
        /// </summary>
        public void IndicateTaskStart()
        {
            _synchronizationContext?.Post(new SendOrPostCallback(state => TaskStarted?.Invoke()), null);
        }

        /// <summary>
        /// Indicate that the task has ended.
        /// </summary>
        public void IndicateTaskEnded()
        {
            _synchronizationContext?.Post(new SendOrPostCallback(state => TaskEnded?.Invoke()), null);
        }

        /// <summary>
        /// Performs action when progress is reported.
        /// </summary>
        /// <param name="prog">The progress level.</param>
        protected virtual void OnProgressReported(double prog)
        {
        }

        /// <summary>
        /// Performs action when message is reported.
        /// </summary>
        /// <param name="msg">The message.</param>
        protected virtual void OnMessageReported(MessageContentStruct msg)
        {
        }

        /// <summary>
        /// Report external process that is running.
        /// </summary>
        /// <param name="process">The external process.</param>
        public void ReportExternalProcess(Process process)
        {
            _externalProcess = process;
            OnExternalProcessReported(process);
        }

        /// <summary>
        /// Performs action when external process is reported.
        /// </summary>
        /// <param name="process">The external process.</param>
        protected virtual void OnExternalProcessReported(Process process)
        {
        }

        /// <summary>
        /// Report the error message.
        /// </summary>
        /// <param name="message">The error message.</param>
        public void ReportError(string message)
        {
            ReportMessage(message, MessageType.FatalError);
        }

        /// <summary>
        /// Invokes the progress handlers.
        /// </summary>
        /// <param name="state">The object.</param>
        private void InvokeProgressHandlers(object state)
        {
            double prog = ((double[])state)[0];
            double prevProg = ((double[])state)[1];
            if (prevProg < 0d)
                prevProg = 0d;
            OnProgressReported(prog);
            ProgressReported?.Invoke(this, prog, prog - prevProg);
        }

        /// <summary>
        /// Invokes the message handlers.
        /// </summary>
        /// <param name="state">The object.</param>
        private void InvokeMessageHandlers(object state)
        {
            MessageContentStruct prog = (MessageContentStruct)state;
            OnMessageReported(prog);
            MessageReported?.Invoke(prog);
        }

        /// <summary>
        /// Cancel the task.
        /// </summary>
        public void RequestCancel()
        {
            _cancellationTokenSource.Cancel();
        }

        /// <summary>
        /// Returns a new thread-safe progress reporter for a subtask.
        /// </summary>
        /// <param name="fractionOfTotal">The fraction of total progress to track.</param>
        /// <param name="subTaskName">The subtask name.</param>
        public SafeProgressReporter CreateProgressModifier(float fractionOfTotal, string subTaskName = "")
        {
            if (fractionOfTotal > 1f)
                fractionOfTotal = 1f;
            if (fractionOfTotal < 0f)
                fractionOfTotal = 0f;
            if (string.IsNullOrEmpty(subTaskName))
                subTaskName = TaskName;
            var child = new SafeProgressReporter(subTaskName);
            child.SetContext(_synchronizationContext);
            child._previousProgress = 0d;
            child.ProgressReported += (reporter, prog, progDelta) => ReportProgress(_previousProgress + progDelta * fractionOfTotal);
            child.MessageReported += msg => ReportMessage(msg);
            child._cancellationTokenSource = _cancellationTokenSource;
            _subProgReporterCollection.Add(child);
            var invokeChildCreatedHandlers = new SendOrPostCallback(state => ChildReporterCreated?.Invoke(child));
            _synchronizationContext?.Post(invokeChildCreatedHandlers, child);
            return child;
        }

        /// <summary>
        /// Report progress and message.
        /// </summary>
        /// <param name="progress">The progress to report.</param>
        /// <param name="message">The message text.</param>
        /// <param name="messageType">The message type.</param>
        public void Report(double progress, string message, MessageType messageType)
        {
            ReportProgress(progress);
            ReportMessage(message, messageType);
        }

        /// <summary>
        /// Report progress.
        /// </summary>
        /// <param name="progress">The progress to report.</param>
        public void ReportProgress(double progress)
        {
            if (_previousProgress == progress) return;
            _synchronizationContext?.Post(_invokeProgressHandlers, new double[] { progress, _previousProgress });
            _previousProgress = progress;
        }

        /// <summary>
        /// Report message.
        /// </summary>
        /// <param name="message">The message text.</param>
        /// <param name="messageType">The message type.</param>
        public void ReportMessage(string message, MessageType messageType = MessageType.Status)
        {
            _synchronizationContext?.Post(_invokeMessageHandlers, new MessageContentStruct(message, messageType, this));
            MessageCount += 1;
            _previousMessage = message;
            _previousMessageType = messageType;
        }

        /// <summary>
        /// Report message.
        /// </summary>
        /// <param name="message">The message to report.</param>
        protected void ReportMessage(MessageContentStruct message)
        {
            _synchronizationContext?.Post(_invokeMessageHandlers, message);
            MessageCount += 1;
            _previousMessage = message.Message;
            _previousMessageType = message.MessageType;
        }

        #endregion

    }
}
