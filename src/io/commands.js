import { spawn } from 'node:child_process';

export function runCommand(command, args) {
	return new Promise((resolve, reject) => {
		const child = spawn(command, args, { stdio: 'inherit' });
		child.on('error', reject);
		child.on('close', (code) => {
			if (code === 0) resolve();
			else reject(new Error(`${command} beendet mit Exit-Code ${code}`));
		});
	});
}
