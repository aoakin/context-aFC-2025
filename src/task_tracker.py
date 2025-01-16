from datetime import datetime

# dict to store task start and end times

task_times = {}

def startTask(task_name):
	current_time = datetime.now()
	current_ISO_time = current_time.ctime()
	task_times[task_name] = {'start': current_time}
	print(f"Starting {task_name} | {current_ISO_time}")

def endTask(task_name):
	current_time = datetime.now()
	current_ISO_time = current_time.ctime()
	if task_name in task_times:
		task_times[task_name]['end'] = current_time
		print(f"Finished {task_name} | {current_ISO_time}")
	else:
		print(f"Task {task_name} was not started.")

def get_elapsed_time(start_time, end_time):
	if start_time and end_time:
		elapsed = end_time - start_time
		return str(elapsed)
	return "N/A"

def TaskTimes():
	for task_name, times in task_times.items():
		start_time = times.get('start', None)
		end_time = times.get('end', None)
		elapsed_time = get_elapsed_time(start_time, end_time)

		print(f"Task: {task_name}")
		print(f"\tStart Time: {start_time.ctime()}")
		print(f"\tEnd Time: {end_time.ctime()}")
		print(f"\tElapsed Time: {elapsed_time}\n")
