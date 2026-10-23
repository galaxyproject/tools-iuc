import { defineConfig } from 'astro/config';
import tailwindcss from '@tailwindcss/vite';
import sitemap from '@astrojs/sitemap';

// The site is served from GitHub Pages at <org>.github.io/tools-iuc
export default defineConfig({
  site: 'https://galaxyproject.github.io',
  base: '/tools-iuc',
  output: 'static',
  integrations: [sitemap()],
  trailingSlash: 'ignore',
  vite: {
    plugins: [tailwindcss()],
  },
});
