using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Utilities;
using System.Reflection;
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
        public void Test_CreateProgressModifier_ResetHandoffIsAtomic()
        {
            var parent = new SafeProgressReporter("parent");
            var lockField = typeof(SafeProgressReporter).GetField("_subProgReporterLock", BindingFlags.Instance | BindingFlags.NonPublic);
            var sourceField = typeof(SafeProgressReporter).GetField("_cancellationTokenSource", BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.IsNotNull(lockField);
            Assert.IsNotNull(sourceField);

            object registryLock = lockField.GetValue(parent)!;
            SafeProgressReporter child = null!;
            var registrationThread = new Thread(() => child = parent.CreateProgressModifier(1f, "child"));
            bool registrationBlocked;

            Monitor.Enter(registryLock);
            try
            {
                registrationThread.Start();
                registrationBlocked = SpinWait.SpinUntil(
                    () => (registrationThread.ThreadState & ThreadState.WaitSleepJoin) != 0,
                    5000);
                if (registrationBlocked)
                    sourceField.SetValue(parent, new CancellationTokenSource());
            }
            finally
            {
                Monitor.Exit(registryLock);
            }

            Assert.IsTrue(registrationBlocked, "Child registration did not reach the registry lock.");
            Assert.IsTrue(registrationThread.Join(5000), "Child registration did not complete.");

            parent.RequestCancel();
            Assert.IsTrue(child.IsCancelRequested, "The child retained the cancellation source from before the reset handoff.");
        }
    }
}
