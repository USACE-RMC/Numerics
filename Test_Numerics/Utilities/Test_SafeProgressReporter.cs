using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Utilities;
using System.Reflection;
using System.Runtime.ExceptionServices;
using System.Threading;
using System.Threading.Tasks;

namespace Utilities
{
    /// <summary>
    /// Unit tests for the thread safety of the progress reporter's child registry.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_SafeProgressReporter
    {
        /// <summary>
        /// Creating child reporters concurrently from a parallel loop registers every child and
        /// links every child to the parent's cancellation source.
        /// </summary>
        [TestMethod]
        public void Test_CreateProgressModifier_ParallelRegistration()
        {
            var parent = new SafeProgressReporter("parent");
            int children = 1000;
            var created = new SafeProgressReporter[children];

            Parallel.For(0, children, i =>
            {
                created[i] = parent.CreateProgressModifier(1f / children, $"child {i}");
            });

            Assert.HasCount(children, parent.ChildReporters);

            parent.RequestCancel();
            for (int i = 0; i < children; i++)
            {
                Assert.IsTrue(created[i].IsCancelRequested, $"child {i} did not receive the cancellation request");
            }
        }

        /// <summary>
        /// Enumerating the child reporters while new children are being registered neither throws
        /// nor observes a partially updated registry entry.
        /// </summary>
        [TestMethod]
        public void Test_ChildReporters_SnapshotDuringRegistration()
        {
            var parent = new SafeProgressReporter("parent");
            int children = 500;

            var addTask = Task.Run(() =>
            {
                for (int i = 0; i < children; i++)
                    parent.CreateProgressModifier(1f / children, $"child {i}");
            });

            while (!addTask.IsCompleted)
            {
                foreach (var child in parent.ChildReporters)
                {
                    Assert.IsNotNull(child);
                }
            }
            addTask.Wait();
            Assert.HasCount(children, parent.ChildReporters);
        }

        /// <summary>
        /// Registering a child across a cancellation reset must link the child to the replacement
        /// cancellation source rather than the source that was current before registration.
        /// </summary>
        [TestMethod]
        [DoNotParallelize]
        public void Test_CreateProgressModifier_ResetHandoffIsAtomic()
        {
            var parent = new SafeProgressReporter("parent");
            var lockField = typeof(SafeProgressReporter).GetField("_subProgReporterLock", BindingFlags.Instance | BindingFlags.NonPublic);
            var sourceField = typeof(SafeProgressReporter).GetField("_cancellationTokenSource", BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.IsNotNull(lockField);
            Assert.IsNotNull(sourceField);

            object registryLock = lockField.GetValue(parent)!;
            SafeProgressReporter child = null!;
            using var registrationWaiting = new ManualResetEventSlim();
            var registrationContext = new WaitNotifyingSynchronizationContext(registrationWaiting);
            ExceptionDispatchInfo registrationException = null!;
            var registrationThread = new Thread(() =>
            {
                var previousContext = SynchronizationContext.Current;
                try
                {
                    SynchronizationContext.SetSynchronizationContext(registrationContext);
                    child = parent.CreateProgressModifier(1f, "child");
                }
                catch (Exception ex)
                {
                    registrationException = ExceptionDispatchInfo.Capture(ex);
                }
                finally
                {
                    SynchronizationContext.SetSynchronizationContext(previousContext);
                }
            }) { IsBackground = true };
            bool registrationStarted = false;
            bool registrationBlocked = false;
            bool registrationCompleted = false;

            Monitor.Enter(registryLock);
            try
            {
                registrationThread.Start();
                registrationStarted = true;
                // The worker performs no other waits after installing the context, so this
                // notification comes from Monitor.Enter inside CreateProgressModifier.
                registrationBlocked = registrationWaiting.Wait(5000);
                if (registrationBlocked)
                    sourceField.SetValue(parent, new CancellationTokenSource());
            }
            finally
            {
                Monitor.Exit(registryLock);
                if (registrationStarted)
                    registrationCompleted = registrationThread.Join(5000);
            }

            Assert.IsTrue(registrationCompleted, "Child registration did not complete.");
            registrationException?.Throw();
            Assert.IsTrue(registrationBlocked, "Child registration did not reach the registry lock.");

            parent.RequestCancel();
            Assert.IsTrue(child.IsCancelRequested, "The child retained the cancellation source from before the reset handoff.");
        }

        /// <summary>
        /// Signals when the registration thread enters a blocking wait, then preserves the
        /// normal wait behavior so the test can replace the source while holding the registry lock.
        /// </summary>
        private sealed class WaitNotifyingSynchronizationContext : SynchronizationContext
        {
            private readonly ManualResetEventSlim _waiting;

            /// <summary>
            /// Creates a context that notifies the test when a blocking wait begins.
            /// </summary>
            /// <param name="waiting">The event to signal when the worker begins waiting.</param>
            public WaitNotifyingSynchronizationContext(ManualResetEventSlim waiting)
            {
                _waiting = waiting;
                SetWaitNotificationRequired();
            }

            /// <summary>
            /// Signals the blocking wait and delegates to the default synchronization context.
            /// </summary>
            /// <param name="waitHandles">The native handles to wait on.</param>
            /// <param name="waitAll">Whether all handles must be signaled.</param>
            /// <param name="millisecondsTimeout">The maximum wait duration in milliseconds.</param>
            /// <returns>The result of the default synchronization-context wait.</returns>
            public override int Wait(IntPtr[] waitHandles, bool waitAll, int millisecondsTimeout)
            {
                _waiting.Set();
                return base.Wait(waitHandles, waitAll, millisecondsTimeout);
            }
        }
    }
}
