import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './tests',
  use: {
    baseURL: 'http://127.0.0.1:4326/tools-iuc/',
    launchOptions: { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH },
  },
  webServer: {
    command: 'npm run preview -- --host 127.0.0.1 --port 4326 --strictPort',
    url: 'http://127.0.0.1:4326/tools-iuc/',
    timeout: 30_000,
  },
});
